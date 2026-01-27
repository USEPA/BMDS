import math
import re
import warnings
from itertools import chain
from pathlib import Path

import numpy as np
import pandas as pd
from pydantic import BaseModel

from ..constants import (
    ContinuousModel,
    DichotomousModel,
    DistType,
    Dtype,
    NestedDichotomousModel,
    PriorClass,
    PriorDistribution,
)
from ..utils import pretty_table


class Prior(BaseModel):
    name: str
    type: PriorDistribution
    initial_value: float
    stdev: float
    min_value: float
    max_value: float

    def numeric_list(self) -> list[float]:
        return list(self.model_dump(exclude={"name"}).values())


class ModelPriors(BaseModel):
    prior_class: PriorClass  # if this is a predefined model class
    priors: list[Prior]  # priors for main model
    variance_priors: list[Prior] | None = None  # priors for variance model (continuous-only)
    overrides: dict[str, dict] | None = None  # beta term overrides

    def report_tbl(self) -> str:
        """Generate a table of priors given this configuration.

        Note that this doesn't include any beta overrides, expansion of polynomial terms,
        or adjustments based on the continuous variance configurations. To get the priors
        used in a specific model, use the model `priors_tbl` method. This method is primarily
        used to investigate default model configurations, without additional settings applied.
        """
        headers = "name|type|initial|stdev|min|max".split("|")
        rows = [
            (p.name, p.type.name, p.initial_value, p.stdev, p.min_value, p.max_value)
            for p in chain(self.priors, self.variance_priors or ())
        ]
        return pretty_table(rows, headers)

    def apply_continuous_loud_defaults(self, dataset, dist_type: DistType) -> None:
        """
        Apply dataset-informed default priors for Continuous LOUD models.
        This overwrites CSV defaults BEFORE user overrides are applied.
        """
        if self.prior_class is not PriorClass.bayesian_loud:
            return

        doses = np.asarray(dataset.doses, dtype=float)
        ns = np.asarray(dataset.ns, dtype=float)
        means = np.asarray(dataset.means, dtype=float)
        stdevs = np.asarray(dataset.stdevs, dtype=float)

        i0 = int(np.argmin(doses))
        i1 = int(np.argmax(doses))

        n0, n1 = ns[i0], ns[i1]
        m0, m1 = means[i0], means[i1]
        s0, s1 = stdevs[i0], stdevs[i1]

        def _p(name: str) -> Prior | None:
            try:
                return self.get_prior(name)
            except ValueError:
                return None

        # -----------------------------
        # Mean priors (Student-t)
        # -----------------------------
        df0 = n0 - 1
        df1 = n1 - 1

        if dist_type in {DistType.normal, DistType.normal_ncv}:
            scale0 = s0 / math.sqrt(df0)
            scale1 = s1 / math.sqrt(df1)

            if p := _p("m0"):
                p.type = PriorDistribution.Student_t
                p.initial_value = df0
                p.stdev = m0  # loc
                p.min_value = scale0  # scale

            if p := _p("m1"):
                p.type = PriorDistribution.Student_t
                p.initial_value = df1
                p.stdev = m1
                p.min_value = scale1

        elif dist_type is DistType.log_normal:
            if m0 <= 0 or m1 <= 0:
                return

            varlog0 = math.log(1 + (s0**2 / m0**2))
            varlog1 = math.log(1 + (s1**2 / m1**2))
            s01 = ((df0 * varlog0) + (df1 * varlog1)) / (df0 + df1)

            scale0 = math.sqrt(s01 / df0)
            scale1 = math.sqrt(s01 / df1)

            if p := _p("m0"):
                p.type = PriorDistribution.Student_t
                p.initial_value = df0
                p.stdev = math.log(m0)
                p.min_value = scale0

            if p := _p("m1"):
                p.type = PriorDistribution.Student_t
                p.initial_value = df1
                p.stdev = math.log(m1)
                p.min_value = scale1

        # -----------------------------
        # Variance priors (Inv-Gamma)
        # -----------------------------
        if dist_type is DistType.normal:
            # pooled variance
            s2 = ((df0 * s0**2) + (df1 * s1**2)) / (df0 + df1)
            if p := _p("Var0"):
                p.type = PriorDistribution.InverseGamma
                p.initial_value = (df0 + df1) / 2  # shape
                p.stdev = s2 * (df0 + df1) / 2  # scale

        elif dist_type is DistType.normal_ncv:
            if p := _p("Var0"):
                p.type = PriorDistribution.InverseGamma
                p.initial_value = df0 / 2
                p.stdev = (s0**2) * df0 / 2

            if p := _p("Var1"):
                p.type = PriorDistribution.InverseGamma
                p.initial_value = df1 / 2
                p.stdev = (s1**2) * df1 / 2

        elif dist_type is DistType.log_normal:
            s2log = ((df0 * varlog0) + (df1 * varlog1)) / (df0 + df1)
            if p := _p("Var0"):
                p.type = PriorDistribution.InverseGamma
                p.initial_value = (df0 + df1) / 2
                p.stdev = s2log * (df0 + df1) / 2

    def get_prior(self, name: str) -> Prior:
        """Search all priors and return the match by name.

        Args:
            name (str): prior name

        Raises:
            ValueError: if no value is found
        """
        for p in chain(self.priors, self.variance_priors or []):
            if p.name == name:
                return p
        raise ValueError(f"No parameter named {name}")

    def update(self, name: str, **kw):
        """Update a prior inplace.

        Args:
            name (str): the prior name
            **kw: fields to update
        """

        # If the term being adjusted is a beta term from a polynomial model; save in the beta
        # overrides instead of altering directly (the polynomial prior expansion is a special case)
        match = re.search(r"^b[2-9]$", name)
        if match:
            if self.overrides is None:
                self.overrides = {}
            self.overrides[match[0]] = kw
            return

        # If the term being adjusted is a phi term from a nested dichotomous model; save in the beta
        # overrides instead of altering directly (the polynomial prior expansion is a special case)
        match = re.search(r"^phi[1-9]$", name)
        if match:
            if self.overrides is None:
                self.overrides = {}
            self.overrides[match[0]] = kw
            return

        # otherwise set revisions directly
        prior = self.get_prior(name)
        if prior.type is PriorDistribution.Student_t:
            alias = {"df": "initial_value", "loc": "stdev", "scale": "min_value"}
        elif prior.type is PriorDistribution.InverseGamma:
            alias = {"shape": "initial_value", "scale": "stdev"}
        else:
            alias = {}

        kw = {alias.get(k, k): v for k, v in kw.items()}

        for k, v in kw.items():
            setattr(prior, k, v)

    def priors_list(
        self,
        degree: int | None = None,
        dist_type: DistType | None = None,
        nphi: int | None = None,
    ) -> list[list]:
        priors = []
        for prior in self.priors:
            if nphi is not None and prior.name == "phi":
                continue
            priors.append(prior.model_copy())

        if degree:
            priors.pop(2)

        # copy degree N; > 2nd order poly
        if degree and degree >= 2:
            overrides = self.overrides or {}
            for i in range(2, degree + 1):
                prior = self.priors[2].model_copy()
                for key, value in overrides.get(f"b{i}", {}).items():
                    setattr(prior, key, value)
                priors.append(prior)

        # copy phi N times
        if nphi:
            overrides = self.overrides or {}
            phi = self.get_prior("phi")
            for i in range(1, nphi + 1):
                prior = phi.model_copy()
                for key, value in overrides.get(f"phi{i}", {}).items():
                    setattr(prior, key, value)
                priors.append(prior)

        if self.prior_class is PriorClass.bayesian_loud:
            if dist_type in {DistType.normal, DistType.log_normal}:
                priors.append(self.get_prior("Var0").model_copy())
            elif dist_type is DistType.normal_ncv:
                priors.append(self.get_prior("Var0").model_copy())
                priors.append(self.get_prior("Var1").model_copy())

        else:
            # add constant variance parameter
            if dist_type and dist_type in {DistType.normal, DistType.log_normal}:
                priors.append(self.variance_priors[1].model_copy())

            # add non-constant variance parameter
            if dist_type and dist_type is DistType.normal_ncv:
                for variance_prior in self.variance_priors:
                    priors.append(variance_prior)

        # check values
        for prior in priors:
            if prior.min_value > prior.max_value:
                warnings.warn(f"Min Value > Max Value ({prior})", stacklevel=2)
            elif prior.initial_value < prior.min_value:
                warnings.warn(f"Initial Value < Min Value ({prior})", stacklevel=2)
            elif prior.initial_value > prior.max_value:
                warnings.warn(f"Initial Value > Max Value ({prior})", stacklevel=2)

        return [prior.numeric_list() for prior in priors]

    def to_c(self, degree: int | None = None, dist_type: DistType | None = None) -> np.ndarray:
        priors = self.priors_list(degree, dist_type)
        return np.array(priors, dtype=np.float64).flatten("F")

    def to_c_nd(self, n_phi: int) -> np.ndarray:
        # Nested dichotomous output C struct only has two columns instead of all 5
        priors = self.priors_list(nphi=n_phi)
        return np.array(priors, dtype=np.float64)[:, 3:].flatten("F")

    @property
    def is_bayesian(self) -> bool:
        return self.prior_class.is_bayesian


# lazy mapping; saves copy as requested
_model_priors: dict[str, ModelPriors] = {}


def _load_model_priors():
    # lazy load model priors from CSV file
    def set_param_type(df):
        df = df.assign(variance_param=False)
        legacy_variance = {"rho", "alpha", "log-alpha"}
        loud_variance = {"Var0", "Var1"}
        df.loc[(df.data_class == "C") & (df.name.isin(legacy_variance)), "variance_param"] = True

        loud_value = PriorClass.bayesian_loud.value
        df.loc[
            (df.data_class == "C") & (df.prior_class == loud_value) & (df.name.isin(loud_variance)),
            "variance_param",
        ] = True
        return df

    def build_priors(df):
        priors = {}
        for (data_class, model_id, prior_class), params in df:
            key = f"{data_class}-{model_id}-{prior_class}"
            gof_priors = params[params.variance_param == False]  # noqa: E712
            var_priors = params[params.variance_param == True]  # noqa: E712
            priors[key] = ModelPriors(
                prior_class=prior_class,
                priors=gof_priors.to_dict("records"),
                variance_priors=var_priors.to_dict("records") if var_priors.shape[0] > 0 else None,
            )
        return priors

    filename = Path(__file__).parent / "priors.csv"
    priors = (
        pd.read_csv(str(filename))
        .pipe(set_param_type)
        .groupby(["data_class", "model_id", "prior_class"])
        .pipe(build_priors)
    )
    _model_priors.update(priors)


def get_dichotomous_prior(model: DichotomousModel, prior_class: PriorClass) -> ModelPriors:
    if len(_model_priors) == 0:
        _load_model_priors()
    key = f"{Dtype.DICHOTOMOUS.value}-{model.id}-{prior_class}"
    return _model_priors[key].model_copy(deep=True)


def get_continuous_prior(
    model: ContinuousModel, prior_class: PriorClass, dataset=None, dist_type: DistType | None = None
) -> ModelPriors:
    if len(_model_priors) == 0:
        _load_model_priors()
    key = f"{Dtype.CONTINUOUS.value}-{model.id}-{prior_class}"
    mp = _model_priors[key].model_copy(deep=True)

    if dataset is not None and prior_class is PriorClass.bayesian_loud:
        mp.apply_continuous_loud_defaults(dataset, dist_type)

    return mp


def get_nested_dichotomous_prior(
    model: NestedDichotomousModel, prior_class: PriorClass
) -> ModelPriors:
    if len(_model_priors) == 0:
        _load_model_priors()
    key = f"{Dtype.NESTED_DICHOTOMOUS.value}-{model.id}-{prior_class}"
    return _model_priors[key].model_copy(deep=True)


def priors_tbl(
    params: list[str],
    priors: list[list],
    is_bayesian: bool,
    dist_type: DistType | None = None,
) -> str:
    """NOTE: values is [type, initial_value, stdev, min_value, max_value]
    For LOUD we interpret these fields as distribution parameters:
    Student_t: (df, loc, scale, min)  => (initial, stdev, min, max)
    InvGamma: (shape, scale, min, max) => (initial, stdev, min, max)
    """
    headers = []
    rows = []
    if is_bayesian:
        headers = "Parameter|Distribution|Definition"
        for name, values in zip(params, priors, strict=True):
            dist = values[0]

            if dist.name.lower() in {"student_t"}:
                df_ = values[1]
                loc = values[2]
                scale = values[3]
                definition = f"df={df_}, loc={loc}, scale={scale}, min={values[4]}"
            elif dist.name.lower() in {"inversegamma"}:
                shape = values[1]
                scale = values[2]
                definition = f"shape={shape}, scale={scale}, min={values[3]}, max={values[4]}"
            else:
                definition = (
                    f"initial={values[1]}, stdev={values[2]}, min={values[3]}, max={values[4]}"
                )
            display_name = name
            if name == "Var0" and dist_type == DistType.normal:
                display_name = "Var"
            elif name == "Var0" and dist_type == DistType.log_normal:
                display_name = "Var_log"
            rows.append((display_name, dist.name, definition))
    else:
        headers = "Parameter|Initial|Min|Max"
        for name, values in zip(params, priors, strict=True):
            rows.append((name, values[1], values[3], values[4]))
    return pretty_table(rows, headers.split("|"))


def multistage_cancer_prior() -> ModelPriors:
    # fmt: off
    priors = [
        Prior(name="g",  type=PriorDistribution.Uniform, initial_value=-17, stdev=0, min_value=-18, max_value=18),
        Prior(name="b1", type=PriorDistribution.Uniform, initial_value=0.1, stdev=0, min_value=0, max_value=1e4),
        Prior(name="bN", type=PriorDistribution.Uniform, initial_value=0.1, stdev=0, min_value=0, max_value=1e4),
    ]
    # fmt: on
    return ModelPriors(
        prior_class=PriorClass.frequentist_restricted, priors=priors, variance_priors=None
    )
