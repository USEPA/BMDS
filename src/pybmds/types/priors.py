import re
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

        # add constant variance parameter
        if dist_type and dist_type in {DistType.normal, DistType.log_normal}:
            priors.append(self.variance_priors[1].model_copy())

        # add non-constant variance parameter
        if dist_type and dist_type is DistType.normal_ncv:
            for variance_prior in self.variance_priors:
                priors.append(variance_prior)
        return [prior.numeric_list() for prior in priors]

    def to_c(self, degree: int | None = None, dist_type: DistType | None = None) -> np.ndarray:
        priors = self.priors_list(degree, dist_type)
        return np.array(priors, dtype=np.float64).flatten("F")

    @property
    def is_bayesian(self) -> bool:
        return self.prior_class.is_bayesian


# lazy mapping; saves copy as requested
_model_priors: dict[str, ModelPriors] = {}


def _load_model_priors():
    # lazy load model priors from CSV file
    def set_param_type(df):
        df = df.assign(variance_param=False)
        df.loc[
            (df.data_class == "C") & (df.name.isin(["rho", "alpha", "log-alpha"])), "variance_param"
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


def get_continuous_prior(model: ContinuousModel, prior_class: PriorClass) -> ModelPriors:
    if len(_model_priors) == 0:
        _load_model_priors()
    key = f"{Dtype.CONTINUOUS.value}-{model.id}-{prior_class}"
    return _model_priors[key].model_copy(deep=True)


def get_nested_dichotomous_prior(
    model: NestedDichotomousModel, prior_class: PriorClass
) -> ModelPriors:
    if len(_model_priors) == 0:
        _load_model_priors()
    key = f"{Dtype.NESTED_DICHOTOMOUS.value}-{model.id}-{prior_class}"
    return _model_priors[key].model_copy(deep=True)


def priors_tbl(params: list[str], priors: list[list], is_bayesian: bool) -> str:
    headers = []
    rows = []
    if is_bayesian:
        headers = "Parameter|Distribution|Initial|Stdev|Min|Max"
        for name, values in zip(params, priors, strict=True):
            rows.append((name, values[0].name, values[1], values[2], values[3], values[4]))
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
