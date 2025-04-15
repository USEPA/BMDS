import warnings
from enum import IntEnum
from random import randrange
from typing import NamedTuple, Self

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from .. import bmdscore, constants
from ..datasets import NestedDichotomousDataset
from ..utils import camel_to_title, multi_lstrip, pretty_table, unique_items
from .common import BOUND_FOOTNOTE, NumpyFloatArray, clean_array, inspect_cpp_obj
from .priors import ModelPriors, PriorClass


class RiskType(IntEnum):
    AddedRisk = 0
    ExtraRisk = 1


class LitterSpecificCovariate(IntEnum):
    Unused = 0
    OverallMean = 1
    ControlGroupMean = 2

    def label(self, fixed_lsc: float | None) -> str:
        text = camel_to_title(self.name)
        if self == self.Unused:
            return text
        lsc = f" ({fixed_lsc:.3f})" if fixed_lsc else ""
        return f"{text}{lsc}"

    @property
    def text(self) -> str:
        return "lsc-" if self == self.Unused else "lsc+"


class IntralitterCorrelation(IntEnum):
    Zero = 0
    Estimate = 1

    @property
    def label(self) -> str:
        return camel_to_title(self.name)

    @property
    def text(self) -> str:
        return "ilc+" if self == self.Estimate else "ilc-"


_bmr_text_map = {
    RiskType.ExtraRisk: "{:.0%} Extra Risk",
    RiskType.AddedRisk: "{:.0%} Added Risk",
}


class NestedDichotomousModelSettings(BaseModel):
    bmr_type: RiskType = RiskType.ExtraRisk
    bmr: float = Field(default=0.1, gt=0)
    alpha: float = Field(default=0.05, gt=0, lt=1)
    litter_specific_covariate: LitterSpecificCovariate = LitterSpecificCovariate.OverallMean
    intralitter_correlation: IntralitterCorrelation = IntralitterCorrelation.Estimate
    estimate_background: bool = True
    bootstrap_iterations: int = Field(default=1000, ge=10, le=1_000_000)
    bootstrap_seed: int = Field(default_factory=lambda: randrange(0, 1_000), ge=0, le=1_000)  # noqa: S311
    bootstrap_n: int = Field(default=3, ge=1, le=10)
    name: str = ""  # override model name
    priors: PriorClass | ModelPriors | None = None  # if None; default used

    model_config = ConfigDict(extra="forbid")

    @property
    def bmr_text(self) -> str:
        return _bmr_text_map[self.bmr_type].format(self.bmr)

    @property
    def confidence_level(self) -> float:
        return 1.0 - self.alpha

    @property
    def modeling_approach(self) -> str:
        return "MLE"

    def tbl(self, results=None) -> str:
        fixed_lsc = results.fixed_lsc if results else None
        return pretty_table(
            [
                ["BMR", self.bmr_text],
                ["Confidence Level (one sided)", self.confidence_level],
                ["Litter Specific Covariate", self.litter_specific_covariate.label(fixed_lsc)],
                ["Intralitter Correlation", self.intralitter_correlation.label],
                ["Estimate Background", self.estimate_background],
                ["Bootstrap Runs", self.bootstrap_n],
                ["Bootstrap Iterations", self.bootstrap_iterations],
                ["Bootstrap Seed", self.bootstrap_seed],
            ],
            "",
        )

    @classmethod
    def docx_table_data(cls, settings: list[Self], results) -> dict:
        fixed_lsc = results.fixed_lsc if results else None
        lsc = [
            setting.litter_specific_covariate.label(fixed_lsc)
            for setting in settings
            if setting.litter_specific_covariate != LitterSpecificCovariate.Unused
        ]
        ilc = [
            setting.intralitter_correlation.label
            for setting in settings
            if setting.intralitter_correlation != IntralitterCorrelation.Zero
        ]
        data = {
            "Setting": "Value",
            "BMR": unique_items(settings, "bmr_text"),
            "Confidence Level (one sided)": unique_items(settings, "confidence_level"),
            "Litter Specific Covariate": ", ".join(sorted(set(lsc))),
            "Intralitter Correlation": ", ".join(sorted(set(ilc))),
            "Estimate Background": unique_items(settings, "estimate_background"),
            "Bootstrap Runs": unique_items(settings, "bootstrap_n"),
            "Bootstrap Seed": unique_items(settings, "bootstrap_seed"),
            "Bootstrap Iterations": unique_items(settings, "bootstrap_iterations"),
        }
        return data

    def update_record(self, d: dict) -> None:
        """Update data record for a tabular-friendly export"""
        d.update(
            bmr=self.bmr_text,
        )


class NestedDichotomousAnalysis(NamedTuple):
    analysis: bmdscore.python_nested_analysis
    result: bmdscore.python_nested_result

    @classmethod
    def create(
        cls,
        model_class: bmdscore.nested_model,
        nparms: int,
        dataset: NestedDichotomousDataset,
        settings: NestedDichotomousModelSettings,
    ) -> Self:
        analysis = bmdscore.python_nested_analysis()
        analysis.model = model_class
        analysis.doses = dataset.doses
        analysis.litterSize = dataset.litter_ns
        analysis.incidence = dataset.incidences
        analysis.lsc = dataset.litter_covariates
        analysis.prior = settings.priors.to_c_nd(n_phi=dataset.num_dose_groups)
        analysis.LSC_type = settings.litter_specific_covariate.value
        analysis.ILC_type = settings.intralitter_correlation.value
        analysis.BMD_type = settings.bmr_type.value
        analysis.estBackground = settings.estimate_background
        analysis.parms = nparms
        analysis.prior_cols = 2
        analysis.BMR = settings.bmr
        analysis.alpha = settings.alpha
        analysis.numBootRuns = settings.bootstrap_n
        analysis.iterations = settings.bootstrap_iterations
        analysis.seed = settings.bootstrap_seed

        result = bmdscore.python_nested_result()
        result.nparms = analysis.parms
        result.model = analysis.model
        result.bmdsRes = bmdscore.BMDS_results()
        result.boot = bmdscore.nestedBootstrap()
        result.litter = bmdscore.nestedLitterData()
        result.reduced = bmdscore.nestedReducedData()
        result.srData = bmdscore.nestedSRData()

        return cls(analysis, result)

    def execute(self):
        bmdscore.pythonBMDSNested(self.analysis, self.result)

    def __str__(self) -> str:
        lines = []
        inspect_cpp_obj(lines, self.analysis, depth=0)
        inspect_cpp_obj(lines, self.result, depth=0)
        return "\n".join(lines)


class BootstrapRuns(BaseModel):
    p_value: list[float]
    p50: list[float]
    p90: list[float]
    p95: list[float]
    p99: list[float]

    @classmethod
    def from_model(cls, data: bmdscore.nestedBootstrap) -> Self:
        return cls(
            p_value=data.pVal,
            p50=data.perc50,
            p90=data.perc90,
            p95=data.perc95,
            p99=data.perc99,
        )

    def tbl(self) -> str:
        col1 = "1 2 3 Combined".split()
        data = list(zip(col1, self.p_value, self.p50, self.p90, self.p95, self.p99, strict=True))
        return pretty_table(data, headers="Run P-Value 50th 90th 95th 99th".split())


class ReducedResult(BaseModel):
    dose: list[float]
    prop_affected: list[float]
    lower_ci: list[float]
    upper_ci: list[float]

    @classmethod
    def from_model(cls, data: bmdscore.nestedReducedData) -> Self:
        return cls(
            dose=data.dose,
            prop_affected=data.propAffect,
            lower_ci=data.lowerConf,
            upper_ci=data.upperConf,
        )


class LitterResult(BaseModel):
    lsc: list[float]
    scaled_residuals: list[float]
    dose: list[float]
    estimated_probabilities: list[float]
    expected: list[float]
    litter_size: list[float]
    observed: list[int]

    def mean_abs_control_residual(self) -> float:
        arr = np.array([self.dose, self.scaled_residuals])
        slice = arr[0] == arr[0].min()
        return np.abs(arr[1, slice]).mean()

    @classmethod
    def from_model(cls, data: bmdscore.nestedLitterData, bmd: float) -> Self:
        return cls(
            lsc=data.LSC,
            scaled_residuals=data.SR,
            dose=data.dose,
            estimated_probabilities=data.estProb,
            expected=data.expected,
            litter_size=data.litterSize,
            observed=data.observed,
        )

    def tbl(self) -> str:
        headers = "Dose|LSC|Est. Prob.|Litter N|Expected|Observed|Scaled Residual".split("|")
        data = list(
            zip(
                self.dose,
                self.lsc,
                self.estimated_probabilities,
                self.litter_size,
                self.expected,
                self.observed,
                self.scaled_residuals,
                strict=True,
            )
        )
        return pretty_table(data, headers)


class Plotting(BaseModel):
    dr_x: NumpyFloatArray
    dr_y: NumpyFloatArray
    bmdl_y: float
    bmd_y: float
    bmdu_y: float

    @classmethod
    def from_model(cls, model, params: dict, fixed_lsc: float) -> Self:
        summary = model.structs.result.bmdsRes
        xs = np.array([summary.BMDL, summary.BMD, summary.BMDU])
        extra_values = [summary.BMD] if summary.BMD >= 0 else []
        dr_x = model.dataset.dose_linspace(extra_values=extra_values)
        dr_y = clean_array(model.dr_curve(dr_x, params, fixed_lsc))
        with warnings.catch_warnings():
            # xs may have a nan if a BMDU is not calculated; filter power warnings
            warnings.filterwarnings("ignore", message=".*invalid value encountered in power.*")
            critical_ys = clean_array(model.dr_curve(xs, params, fixed_lsc))
        critical_ys[critical_ys <= 0] = constants.BMDS_BLANK_VALUE
        return cls(
            dr_x=dr_x,
            dr_y=dr_y,
            bmdl_y=critical_ys[0],
            bmd_y=critical_ys[1],
            bmdu_y=critical_ys[2],
        )


class ScaledResidual(BaseModel):
    min: float
    avg: float
    max: float
    min_abs: float
    avg_abs: float
    max_abs: float

    @classmethod
    def from_model(cls, data: bmdscore.nestedSRData) -> Self:
        return cls(
            min=data.minSR,
            avg=data.avgSR,
            max=data.maxSR,
            min_abs=data.minAbsSR,
            avg_abs=data.avgAbsSR,
            max_abs=data.maxAbsSR,
        )

    def tbl(self) -> str:
        data = [
            ["Minimum scaled residual", self.min],
            ["Minimum ABS(scaled residual)", self.min_abs],
            ["Average scaled residual", self.avg],
            ["Average ABS(scaled residual)", self.avg_abs],
            ["Maximum scaled residual", self.max],
            ["Maximum ABS(scaled residual)", self.max_abs],
        ]
        return pretty_table(data, "")


class NestedDichotomousResult(BaseModel):
    bmd: float
    bmdl: float
    bmdu: float
    aic: float
    bic_equiv: float
    chi_squared: float
    bounded: list[bool]
    lower_ci: list[float]
    std_err: list[float]
    upper_ci: list[float]
    scaled_residuals: ScaledResidual
    bootstrap: BootstrapRuns
    combined_pvalue: float
    ll: float
    cov: list[float]
    dof: float
    fixed_lsc: float
    litter: LitterResult
    parameter_names: list[str]
    parameters: list[float]
    reduced: ReducedResult
    plotting: Plotting
    has_completed: bool = False

    @classmethod
    def from_model(cls, model) -> Self:
        result: bmdscore.python_nested_result = model.structs.result
        params_d = {
            name: value for name, value in zip(model.get_param_names(), result.parms, strict=True)
        }
        return cls(
            bmd=result.bmdsRes.BMD,
            bmdl=result.bmdsRes.BMDL,
            bmdu=result.bmdsRes.BMDU,
            aic=result.bmdsRes.AIC,
            bic_equiv=result.bmdsRes.BIC_equiv,
            bounded=result.bmdsRes.bounded,
            chi_squared=result.bmdsRes.chisq,
            lower_ci=result.bmdsRes.lowerConf,
            std_err=result.bmdsRes.stdErr,
            upper_ci=result.bmdsRes.upperConf,
            scaled_residuals=ScaledResidual.from_model(result.srData),
            bootstrap=BootstrapRuns.from_model(result.boot),
            combined_pvalue=result.combPVal,
            ll=result.LL,
            cov=result.cov,
            dof=result.model_df,
            fixed_lsc=result.fixedLSC,
            litter=LitterResult.from_model(result.litter, result.bmdsRes.BMD),
            parameter_names=list(params_d.keys()),
            parameters=list(params_d.values()),
            reduced=ReducedResult.from_model(result.reduced),
            plotting=Plotting.from_model(model, params_d, result.fixedLSC),
            has_completed=result.validResult,
        )

    def text(
        self, dataset: NestedDichotomousDataset, settings: NestedDichotomousModelSettings
    ) -> str:
        return multi_lstrip(
            f"""
        Modeling Summary:
        {self.tbl()}

        Model Parameters:
        {self.parameter_tbl()}

        Bootstrap Runs:
        {self.bootstrap.tbl()}

        Scaled Residuals (for dose group nearest the BMD):
        {self.scaled_residuals.tbl()}

        Litter Data:
        {self.litter.tbl()}
        """
        )

    def tbl(self) -> str:
        data = [
            ["BMD", self.bmd],
            ["BMDL", self.bmdl],
            ["BMDU", self.bmdu],
            ["AIC", self.aic],
            ["P-Value", self.combined_pvalue],
            ["d.f.", self.dof],
            ["Chi²", self.chi_squared],
            ["Log-Likelihood", self.ll],
        ]
        return pretty_table(data, "")

    def parameter_tbl(self) -> str:
        headers = "Variable|Estimate|On Bound|Std Error".split("|")
        data = []
        for name, value, bounded, se in zip(
            self.parameter_names,
            self.parameters,
            self.bounded,
            self.std_err,
            strict=True,
        ):
            data.append(
                (
                    name,
                    value,
                    constants.BOOL_YES_NO[bounded],
                    "Not Reported" if bounded or se == constants.BMDS_BLANK_VALUE else f"{se:g}",
                )
            )
        text = pretty_table(data, headers)
        if any(self.bounded):
            text += BOUND_FOOTNOTE
        return text

    def parameter_rows(self, extras: dict) -> list[dict]:
        rows = []
        for i in range(len(self.parameter_names)):
            rows.append(
                {
                    **extras,
                    **dict(
                        name=self.parameter_names[i],
                        value=self.parameters[i],
                    ),
                }
            )
        return rows

    def update_record(self, d: dict) -> None:
        """Update data record for a tabular-friendly export"""
        d.update(
            bmd=self.bmd,
            bmdl=self.bmdl,
            bmdu=self.bmdu,
            aic=self.aic,
            p_value=self.combined_pvalue,
            loglikelihood=self.ll,
            dof=self.dof,
            fixed_lsc=self.fixed_lsc,
        )

    def get_parameter(self, parameter: str) -> float:
        """Get parameter value by name"""
        match parameter:
            case "bmd":
                return self.bmd
            case "bmdl":
                return self.bmdl
            case "bmdu":
                return self.bmdu
            case "aic":
                return self.aic
            case "dof":
                return self.dof
            case "pvalue":
                return self.combined_pvalue
            case "roi":
                return self.scaled_residuals.avg_abs
            case "roi_control":
                return self.litter.mean_abs_control_residual()
            case "n_params":
                return len(self.parameters)
            case _:  # pragma: no cover
                raise ValueError(f"Unknown parameter: {parameter}")
