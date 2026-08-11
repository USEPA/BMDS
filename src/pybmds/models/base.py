from __future__ import annotations

import abc
import logging
from typing import TYPE_CHECKING, NamedTuple, Self

import numpy as np
import pandas as pd
from matplotlib.figure import Figure
from pydantic import BaseModel

from .. import plotting
from ..constants import BmdModelSchema as BmdModelClass
from ..constants import ContinuousModelChoices, DistType, Dtype, PriorClass
from ..datasets.base import DatasetType
from ..types.priors import ModelPriors
from ..types.priors import priors_tbl as priors_tbl_fn
from ..utils import get_version, multi_lstrip

if TYPE_CHECKING:
    from ..session import Session


logger = logging.getLogger(__name__)


InputModelSettings = dict | BaseModel | None


def cdf_df(arr: np.ndarray, n_points: int = 200) -> pd.DataFrame:
    arr = np.asarray(arr, dtype=float)

    # Existing BMDS-style CDF format:
    # row 0 = BMD values, row 1 = cumulative probabilities
    if arr.ndim == 2 and 2 in arr.shape:
        values = arr.T if arr.shape[0] == 2 else arr
        percentiles = values[:, 1]
        is_bmds_cdf = (
            np.isfinite(percentiles).all()
            and np.all((0 <= percentiles) & (percentiles <= 1))
            and np.all(np.diff(percentiles) >= 0)
        )
        if is_bmds_cdf:
            df = pd.DataFrame(data=values, columns=["BMD", "Percentile"])
            df["Percentile"] = df["Percentile"] * 100
            return df[["Percentile", "BMD"]]

    # LOUD posterior draws: convert raw draws into an empirical CDF with n_points rows
    # Multiple MCMC chains arrive as chains x iterations and are combined here.
    if arr.ndim in (1, 2):
        draws = arr[np.isfinite(arr)]
        if draws.size == 0:
            return pd.DataFrame(columns=["Percentile", "BMD"])

        percentiles = np.linspace(0.5, 99.5, n_points)
        bmd_values = np.percentile(draws, percentiles)

        return pd.DataFrame(
            {
                "Percentile": percentiles,
                "BMD": bmd_values,
            }
        )

    raise ValueError(f"Unsupported CDF input shape: {arr.shape}")


def cdf_plot(
    cdf: pd.DataFrame,
    alpha: float,
    bmd: float,
    bmdl: float | None,
    bmdu: float | None,
    xlabel: str,
    figsize: tuple[float, float] | None = None,
) -> Figure:
    fig = plotting.create_empty_figure(figsize=figsize)
    ax = fig.gca()
    ax.set_title("BMD Cumulative Distribution Function (CDF)")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Cumulative Probability")
    ax.plot(cdf.BMD, cdf.Percentile / 100, label="Model Probability", **plotting.LINE_FORMAT)
    if bmdl or bmdu:
        if bmdl:
            label = "BMDL, BMDU" if bmdu else "BMDL"
            ax.plot([bmdl, bmdl], [0, alpha], label=label, **plotting.CDF_DASHED)
        if bmdu:
            label = None if bmdl else "BMDU"
            ax.plot([bmdu, bmdu], [0, 1 - alpha], label=label, **plotting.CDF_DASHED)
    if bmd:
        ax.plot([bmd, bmd], [0, 0.5], label="BMD", **plotting.CDF_SOLID)
    ax.legend(**plotting.LEGEND_OPTS)
    return fig


class BmdModel(abc.ABC):
    """
    Captures modeling configuration for model execution.
    Should save no results form model execution or any dataset-specific settings.
    """

    bmd_model_class: BmdModelClass
    degree_required: bool = False

    def __init__(self, dataset: DatasetType, settings: InputModelSettings = None):
        self.dataset = dataset
        self.settings = self.get_model_settings(dataset, settings)
        self.structs: NamedTuple | None = None  # used for model averaging
        self.results: BaseModel | None = None
        self.session: Session | None = None

    def name(self) -> str:
        # return name of model; may be setting-specific
        return self.settings.name or self.bmd_model_class.verbose

    @property
    def has_results(self) -> bool:
        return self.results is not None and self.results.has_completed is True

    @abc.abstractmethod
    def get_model_settings(
        self, dataset: DatasetType, settings: InputModelSettings
    ) -> BaseModel: ...

    @abc.abstractmethod
    def execute(self) -> BaseModel: ...

    def execute_job(self):
        self.execute()

    def _execute_standalone_loud(self) -> BaseModel:
        """Execute a LOUD model through an internal one-model averaging session."""
        from ..session import Session

        original_session = self.session
        session = Session(dataset=self.dataset)
        session.models = [self]
        self.session = session
        try:
            session.add_model_averaging()
            session.execute()
        finally:
            self.session = original_session

        if self.results is None:  # pragma: no cover
            raise RuntimeError("Standalone LOUD execution did not return model results.")
        return self.results

    @abc.abstractmethod
    def serialize(self) -> BaseModel: ...

    def text(self) -> str:
        """Text representation of model inputs and outputs outputs."""
        title = f"{self.name()} Model".center(30) + "\n══════════════════════════════"
        version = get_version()
        version = f"Version: pybmds {version.python} (bmdscore {version.dll})"
        settings = self.model_settings_text()
        if self.has_results:
            results = self.results_text()
        else:
            results = "Model has not successfully executed; no results available."

        return "\n\n".join([title, version, settings, results]) + "\n"

    def results_text(self) -> str:
        if self.results is None:  # pragma: no cover
            raise ValueError("Cannot render text if results are unavailable")
        return self.results.text(self.dataset, self.settings)

    def priors_tbl(self) -> str:
        """Show prior or parameter boundary testing."""
        mp = self.settings.priors
        priors = self.get_priors_list()

        # If priors are LOUD, build the table names from the priors object
        if isinstance(mp, ModelPriors) and mp.prior_class is PriorClass.bayesian_loud:
            names = [p.name for p in mp.priors]
            if getattr(self.bmd_model_class, "id", None) == ContinuousModelChoices.exp_m3.value.id:
                names = [n for n in names if n != "c"]

            var_names = [p.name for p in (mp.variance_priors or [])]
            if getattr(self.settings, "disttype", None) != DistType.normal_ncv and var_names:
                var_names = [var_names[0]]
            names = names + var_names

            return priors_tbl_fn(
                names, priors, mp.is_bayesian, getattr(self.settings, "disttype", None)
            )

        # Default behavior (frequentist or non-LOUD bayesian): use model param names
        return priors_tbl_fn(
            self.get_param_names(), priors, mp.is_bayesian, getattr(self.settings, "disttype", None)
        )

    def model_settings_text(self) -> str:
        return multi_lstrip(
            f"""
        Input Summary:
        {self.settings.tbl(self.degree_required)}

        Parameter Settings:
        {self.priors_tbl()}
        """
        )

    def _plot_bmr_lines(self, ax, axlines: bool):
        plotting.add_bmr_lines(
            ax,
            self.results.bmd,
            self.results.plotting.bmd_y,
            self.results.bmdl,
            self.results.bmdu,
            axlines=axlines,
        )

    def plot(self, figsize: tuple[float, float] | None = None, axlines: bool = False):
        """After model execution, print the dataset, curve-fit, BMD, and BMDL.

        Args:
            figsize (tuple[float, float], optional): Specify an alternative figure size (w, h).
            axlines (bool, optional): Draw the axlines for BMD and BMDL, like legacy BMD plots.
                By default the diamond-based line w/ BMDL and BMDU are used instead.
        """
        if not self.has_results:
            raise ValueError("Cannot plot if results are unavailable")

        fig = self.dataset.plot(figsize=figsize)
        ax = fig.gca()
        if self.dataset.dtype == Dtype.DICHOTOMOUS:
            ax.set_ylim(-0.05, 1.05)
        title = f"{self.dataset._get_dataset_name()}\n{self.name()} Model ({self.settings.modeling_approach})\n{self.settings.bmr_text}"
        ax.set_title(title)
        ax.plot(
            self.results.plotting.dr_x,
            self.results.plotting.dr_y,
            **plotting.LINE_FORMAT,
        )

        # if the BMD is greater than the current plot domain, extend to 105% of BMD
        xlim = ax.get_xlim()
        bmd = self.results.get_parameter("bmd")
        if bmd > xlim[1]:
            ax.set_xlim(xlim[0], bmd * 1.05)

        self._plot_bmr_lines(ax, axlines=axlines)
        slope_factor = getattr(self.results, "slope_factor", None)
        if draw_slope := slope_factor and slope_factor > 0:
            ax.plot(
                [self.results.plotting.dr_x[0], self.results.bmdl],
                [self.results.plotting.dr_y[0], self.results.plotting.bmd_y],
                label="Slope Factor",
                **{**plotting.LINE_FORMAT, "linestyle": (0, (2, 1, 1, 1)), "c": "#ef7215"},
            )

        # reorder handles and labels
        handles, labels = ax.get_legend_handles_labels()
        if draw_slope:
            order = [2, 0, 1] if axlines else [1, 2, 0]
        else:
            order = [1, 0] if axlines else [0, 1]
        ax.legend(
            [handles[idx] for idx in order], [labels[idx] for idx in order], **plotting.LEGEND_OPTS
        )
        plotting.improve_bmd_diamond_legend(ax)
        fig.tight_layout()
        return fig

    def cdf_plot(self, figsize: tuple[float, float] | None = None) -> Figure:
        if not self.has_results:
            raise ValueError("Cannot plot if results are unavailable")
        return cdf_plot(
            cdf=self.cdf(),
            alpha=self.settings.alpha,
            bmd=self.results.bmd,
            bmdl=self.results.bmdl,
            bmdu=self.results.bmdu,
            xlabel=self.dataset.get_xlabel(),
            figsize=figsize,
        )

    def cdf(self) -> pd.DataFrame:
        if not self.has_results:
            raise ValueError("Cannot create if results are unavailable")
        return cdf_df(self.results.fit.bmd_dist)

    @abc.abstractmethod
    def get_param_names(self) -> list[str]: ...

    @abc.abstractmethod
    def get_priors_list(self) -> list[list]: ...

    def to_dict(self) -> dict:
        return self.serialize().model_dump()

    @abc.abstractmethod
    def get_gof_pvalue(self) -> float: ...


class BmdModelSchema(BaseModel):
    @classmethod
    def get_subclass(cls, dtype: Dtype) -> Self:
        from .continuous import BmdModelContinuousSchema
        from .dichotomous import BmdModelDichotomousSchema
        from .nested_dichotomous import BmdModelNestedDichotomousSchema

        if dtype == Dtype.DICHOTOMOUS:
            return BmdModelDichotomousSchema
        elif dtype in Dtype.CONTINUOUS_DTYPES():
            return BmdModelContinuousSchema
        elif dtype == Dtype.NESTED_DICHOTOMOUS:
            return BmdModelNestedDichotomousSchema
        else:
            raise ValueError(f"Invalid dtype: {dtype}")


class BmdModelAveraging(abc.ABC):
    """
    Captures modeling configuration for model execution.
    Should save no results form model execution or any dataset-specific settings.
    """

    def __init__(
        self,
        session: Session,
        models: list[BmdModel],
        settings: InputModelSettings = None,
    ):
        self.session = session
        self.models = models
        # if not settings are not specified copy settings from first model
        initial_settings = settings if settings is not None else models[0].settings
        self.settings = self.get_model_settings(initial_settings)
        self.results: BaseModel | None = None

    @abc.abstractmethod
    def get_model_settings(self, settings: InputModelSettings) -> BaseModel: ...

    @abc.abstractmethod
    def execute(self) -> BaseModel: ...

    def execute_job(self):
        self.results = self.execute()

    @property
    def has_results(self) -> bool:
        return self.results is not None

    @abc.abstractmethod
    def serialize(self, session) -> BmdModelAveragingSchema: ...

    def to_dict(self) -> dict:
        return self.serialize.model_dump()

    def cdf_plot(
        self,
        xlabel: str,
        figsize: tuple[float, float] | None = None,
        cdf: pd.DataFrame | None = None,
    ) -> Figure:
        if not self.has_results:
            raise ValueError("Cannot plot if results are unavailable")
        return cdf_plot(
            cdf=self.cdf() if cdf is None else cdf,
            alpha=self.settings.alpha,
            bmd=self.results.bmd,
            bmdl=self.results.bmdl,
            bmdu=self.results.bmdu,
            xlabel=xlabel,
            figsize=figsize,
        )

    def cdf(self) -> pd.DataFrame:
        if not self.has_results:
            raise ValueError("Cannot create if results are unavailable")
        return cdf_df(self.results.bmd_dist)


class BmdModelAveragingSchema(BaseModel):
    @classmethod
    def get_subclass(cls, dtype: Dtype) -> Self:
        from .cma import BmdModelAveragingContinuousSchema
        from .ma import BmdModelAveragingDichotomousSchema

        if dtype == Dtype.DICHOTOMOUS:
            return BmdModelAveragingDichotomousSchema
        elif dtype in (Dtype.CONTINUOUS, Dtype.CONTINUOUS_INDIVIDUAL):
            return BmdModelAveragingContinuousSchema
        else:
            raise ValueError(f"Invalid dtype: {dtype}")
