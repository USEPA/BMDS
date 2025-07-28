import warnings
from enum import IntEnum
from typing import Annotated, NamedTuple, Self

import numpy as np
from pydantic import BaseModel, ConfigDict, Field

from .. import bmdscore, constants
from ..constants import BOOL_YES_NO, ContinuousModelChoices, Dtype
from ..datasets.continuous import ContinuousDatasets
from ..utils import multi_lstrip, pretty_table, unique_items
from .common import (
    BOUND_FOOTNOTE,
    CONTINUOUS_TEST_FOOTNOTES,
    NumpyFloatArray,
    NumpyIntArray,
    clean_array,
    inspect_cpp_obj,
    residual_of_interest,
)
from .priors import ModelPriors, PriorClass, PriorDistribution


class ContinuousRiskType(IntEnum):
    AbsoluteDeviation = 1
    StandardDeviation = 2
    RelativeDeviation = 3
    PointEstimate = 4
    Extra = 5  # Not used
    HybridExtra = 6
    HybridAdded = 7


_bmr_text_map = {
    ContinuousRiskType.AbsoluteDeviation: "{} Absolute Deviation",
    ContinuousRiskType.StandardDeviation: "{} Standard Deviation",
    ContinuousRiskType.RelativeDeviation: "{:.0%} Relative Deviation",
    ContinuousRiskType.PointEstimate: "{} Point Estimation",
    ContinuousRiskType.Extra: "{} Extra",
    ContinuousRiskType.HybridExtra: "{} Hybrid Extra",
    ContinuousRiskType.HybridAdded: "{} Hybrid Added",
}


class ContinuousModelSettings(BaseModel):
    name: str = ""  # override model name
    bmr_type: ContinuousRiskType = ContinuousRiskType.StandardDeviation
    is_increasing: bool | None = None  # if None; autodetect used
    bmr: float = 1.0
    tail_prob: Annotated[float, Field(gt=0, lt=1)] = 0.01
    disttype: constants.DistType = constants.DistType.normal
    alpha: Annotated[float, Field(gt=0, lt=1)] = 0.05
    samples: Annotated[int, Field(ge=0, le=1000)] = 100
    degree: Annotated[int, Field(ge=0, le=8)] = 0  # polynomial only
    burnin: Annotated[int, Field(ge=5, le=1000)] = 20
    priors: PriorClass | ModelPriors | None = None  # if None; default used
    penalize_aic_on_boundary: bool = True

    model_config = ConfigDict(extra="forbid")

    @property
    def bmr_text(self) -> str:
        return _bmr_text_map[self.bmr_type].format(self.bmr)

    @property
    def direction(self) -> str:
        return "Up (↑)" if self.is_increasing else "Down (↓)"

    @property
    def confidence_level(self) -> float:
        return 1.0 - self.alpha

    @property
    def distribution(self) -> str:
        return f"{self.disttype.distribution_type} + {self.disttype.variance_model}"

    @property
    def is_hybrid(self) -> bool:
        return self.bmr_type in [ContinuousRiskType.HybridExtra, ContinuousRiskType.HybridAdded]

    @property
    def modeling_approach(self) -> str:
        return "Bayesian" if self.priors.is_bayesian else "MLE"

    def tbl(self, show_degree: bool = True) -> str:
        data = [
            ["BMR", self.bmr_text],
            ["Distribution", self.distribution],
            ["Modeling Direction", self.direction],
            ["Confidence Level (one sided)", self.confidence_level],
            ["Modeling Approach", self.modeling_approach],
        ]

        if self.is_hybrid:
            data.append(["Tail Probability", self.tail_prob])

        if show_degree:
            data.append(["Degree", self.degree])

        return pretty_table(data, "")

    @classmethod
    def docx_table_data(cls, settings: list[Self], results) -> dict:
        data = {
            "Setting": "Value",
            "BMR": unique_items(settings, "bmr_text"),
            "Distribution": unique_items(settings, "distribution"),
            "Adverse Direction": unique_items(settings, "direction"),
            "Maximum Polynomial Degree": str(max(setting.degree for setting in settings)),
            "Confidence Level (one sided)": unique_items(settings, "confidence_level"),
        }
        if settings[0].is_hybrid:
            data["Tail Probability"] = unique_items(settings, "tail_prob")
        return data

    def update_record(self, d: dict) -> None:
        """Update data record for a tabular-friendly export"""
        d.update(
            bmr=self.bmr_text,
            distribution=self.distribution,
            direction=self.direction,
            confidence_level=self.confidence_level,
            tail_probability=self.tail_prob,
            degree=self.degree,
            model_class=self.priors.prior_class.name,
        )


MODEL_ENUM_MAP = {
    ContinuousModelChoices.power.value.id: bmdscore.cont_model.power,
    ContinuousModelChoices.hill.value.id: bmdscore.cont_model.hill,
    ContinuousModelChoices.polynomial.value.id: bmdscore.cont_model.polynomial,
    ContinuousModelChoices.exp_m3.value.id: bmdscore.cont_model.exp_3,
    ContinuousModelChoices.exp_m5.value.id: bmdscore.cont_model.exp_5,
}


class ContinuousAnalysis(BaseModel):
    model: constants.ContinuousModel
    dataset: ContinuousDatasets
    priors: ModelPriors
    BMD_type: ContinuousRiskType
    is_increasing: bool
    BMR: float
    tail_prob: float
    disttype: constants.DistType
    alpha: float
    samples: int
    burnin: int
    degree: int
    penalize_aic_on_boundary: bool

    model_config = ConfigDict(arbitrary_types_allowed=True)

    @property
    def num_params(self) -> int:
        if self.model == ContinuousModelChoices.polynomial.value:
            params = self.degree + 1
        else:
            params = len(self.model.params)

        if self.disttype is constants.DistType.normal_ncv:
            params += 2
        else:
            params += 1

        return params

    def _priors_array(self) -> np.ndarray:
        degree = (
            self.degree if self.model.id == ContinuousModelChoices.polynomial.value.id else None
        )
        return self.priors.to_c(degree=degree, dist_type=self.disttype)

    def to_cpp(self):
        analysis = bmdscore.python_continuous_analysis()
        analysis.model = MODEL_ENUM_MAP[self.model.id]
        analysis.BMD_type = self.BMD_type.value
        analysis.BMR = self.BMR
        analysis.parms = self.num_params
        analysis.prior_cols = constants.NUM_PRIOR_COLS
        analysis.transform_dose = 0
        analysis.prior = self._priors_array()
        analysis.degree = self.degree
        analysis.disttype = self.disttype.value
        analysis.alpha = self.alpha
        analysis.tail_prob = self.tail_prob
        analysis.penalizeAIC = self.penalize_aic_on_boundary

        # these 3 variables are related; if setting direction; set others to False
        analysis.isIncreasing = self.is_increasing
        analysis.detectAdvDir = False
        analysis.restricted = False

        if self.dataset.dtype == Dtype.CONTINUOUS:
            analysis.suff_stat = True
            analysis.n = self.dataset.num_dose_groups
            analysis.doses = self.dataset.doses
            analysis.n_group = self.dataset.ns
            analysis.Y = self.dataset.means
            analysis.sd = self.dataset.stdevs
        elif self.dataset.dtype == Dtype.CONTINUOUS_INDIVIDUAL:
            analysis.suff_stat = False
            analysis.n = len(self.dataset.individual_doses)
            analysis.doses = self.dataset.individual_doses
            analysis.n_group = []
            analysis.Y = self.dataset.responses
            analysis.sd = []
        else:
            raise ValueError(f"Invalid dtype: {self.dataset.dtype}")

        result = bmdscore.python_continuous_model_result()
        result.model = analysis.model
        result.dist_numE = constants.N_BMD_DIST
        result.nparms = analysis.parms
        result.gof = bmdscore.continuous_GOF()
        result.bmdsRes = bmdscore.BMDS_results()
        result.aod = bmdscore.continuous_AOD()
        result.aod.TOI = bmdscore.testsOfInterest()

        return ContinuousAnalysisCPPStructs(analysis, result)


class ContinuousAnalysisCPPStructs(NamedTuple):
    analysis: bmdscore.python_continuous_analysis
    result: bmdscore.python_continuous_model_result

    def execute(self):
        bmdscore.pythonBMDSCont(self.analysis, self.result)

    def __str__(self) -> str:
        lines = []
        inspect_cpp_obj(lines, self.analysis, depth=0)
        inspect_cpp_obj(lines, self.result, depth=0)
        return "\n".join(lines)


class ContinuousModelResult(BaseModel):
    dist: int
    loglikelihood: float
    aic: float
    bic_equiv: float
    chisq: float
    bmds_model_df: float = Field(alias="model_df")
    total_df: float
    bmd_dist: NumpyFloatArray

    @classmethod
    def from_model(cls, model) -> Self:
        result = model.structs.result
        summary = result.bmdsRes
        arr = np.array(result.bmd_dist, dtype=float).reshape(2, constants.N_BMD_DIST)
        arr = arr[:, np.isfinite(arr[0, :])]
        arr = arr[:, arr[0, :] > 0]

        return ContinuousModelResult(
            dist=result.dist,
            loglikelihood=result.aod.LL[3],
            aic=summary.AIC,
            bic_equiv=summary.BIC_equiv,
            chisq=summary.chisq,
            model_df=result.model_df,
            total_df=result.total_df,
            bmd_dist=arr,
        )


class ContinuousParameters(BaseModel):
    names: list[str]
    values: NumpyFloatArray
    se: NumpyFloatArray
    lower_ci: NumpyFloatArray
    upper_ci: NumpyFloatArray
    bounded: NumpyFloatArray
    cov: NumpyFloatArray
    prior_type: NumpyIntArray
    prior_initial_value: NumpyFloatArray
    prior_stdev: NumpyFloatArray
    prior_min_value: NumpyFloatArray
    prior_max_value: NumpyFloatArray

    @classmethod
    def get_priors(cls, model) -> np.ndarray:
        priors_list = model.get_priors_list()
        return np.array(priors_list, dtype=np.float64).T

    @classmethod
    def from_model(cls, model) -> Self:
        result = model.structs.result
        summary = result.bmdsRes
        param_names = model.get_param_names()
        priors = cls.get_priors(model)

        cov_n = result.nparms
        cov = np.array(result.cov).reshape(cov_n, cov_n)
        slice = None

        # DLL deletes the c parameter and shifts items down; correct in outputs here
        if model.bmd_model_class.id == ContinuousModelChoices.exp_m3.value.id:
            # do the same for parameter names for consistency
            c_index = param_names.index("c")
            param_names.pop(c_index)

            # shift priors as well
            priors = priors.T
            priors[c_index:-1] = priors[c_index + 1 :]
            priors = priors[:-1].T

            # remove final element for some params (stdErr, lowerConf, upperConf)
            slice = -1

        return cls(
            names=param_names,
            values=result.parms,
            bounded=summary.bounded,
            se=summary.stdErr[:slice],
            lower_ci=summary.lowerConf[:slice],
            upper_ci=summary.upperConf[:slice],
            cov=cov,
            prior_type=priors[0],
            prior_initial_value=priors[1],
            prior_stdev=priors[2],
            prior_min_value=priors[3],
            prior_max_value=priors[4],
        )

    def tbl(self) -> str:
        headers = "Variable|Estimate|On Bound|Std Error".split("|")
        data = []
        for name, value, bounded, se in zip(
            self.names,
            self.values,
            self.bounded,
            self.se,
            strict=True,
        ):
            data.append(
                (
                    name,
                    value,
                    BOOL_YES_NO[bounded],
                    "Not Reported" if bounded else f"{se:g}",
                )
            )
        text = pretty_table(data, headers)
        if any(self.bounded):
            text += BOUND_FOOTNOTE
        return text

    def rows(self, extras: dict) -> list[dict]:
        rows = []
        for i in range(len(self.names)):
            rows.append(
                {
                    **extras,
                    **dict(
                        name=self.names[i],
                        value=self.values[i],
                        se=self.se[i],
                        bounded=bool(self.bounded[i]),
                        initial_distribution=PriorDistribution(self.prior_type[i]).name,
                        initial_value=self.prior_initial_value[i],
                        initial_stdev=self.prior_stdev[i],
                        initial_min_value=self.prior_min_value[i],
                        initial_max_value=self.prior_max_value[i],
                    ),
                }
            )
        return rows


class ContinuousGof(BaseModel):
    dose: NumpyFloatArray
    size: NumpyFloatArray
    est_mean: NumpyFloatArray
    calc_mean: NumpyFloatArray
    obs_mean: NumpyFloatArray
    est_sd: NumpyFloatArray
    calc_sd: NumpyFloatArray
    obs_sd: NumpyFloatArray
    residual: NumpyFloatArray
    eb_lower: NumpyFloatArray
    eb_upper: NumpyFloatArray
    roi: float

    @classmethod
    def from_model(cls, model) -> Self:
        gof = model.structs.result.gof
        summary = model.structs.result.bmdsRes
        return ContinuousGof(
            dose=np.array(gof.dose),
            size=np.array(gof.size),
            est_mean=np.array(gof.estMean),
            calc_mean=np.array(gof.calcMean),
            obs_mean=np.array(gof.obsMean),
            est_sd=np.array(gof.estSD),
            calc_sd=np.array(gof.calcSD),
            obs_sd=np.array(gof.obsSD),
            residual=np.array(gof.res),
            eb_lower=np.array(gof.ebLower),
            eb_upper=np.array(gof.ebUpper),
            roi=residual_of_interest(summary.BMD, model.dataset.doses, np.array(gof.res).tolist()),
        )

    def get_tbl_data_means(self, disttype: constants.DistType):
        if disttype == constants.DistType.log_normal:
            headers = (
                "Dose|N|Sample Mean|Approximate Sample Median|Model Fitted Median|Scaled Residual"
            )
            data = [
                (
                    self.dose[idx],
                    self.size[idx],
                    self.obs_mean[idx],
                    self.calc_mean[idx],
                    self.est_mean[idx],
                    self.residual[idx],
                )
                for idx in range(len(self.dose))
            ]
        else:
            headers = "Dose|N|Sample Mean|Model Fitted Mean|Scaled Residual"
            data = [
                (
                    self.dose[idx],
                    self.size[idx],
                    self.obs_mean[idx],
                    self.est_mean[idx],
                    self.residual[idx],
                )
                for idx in range(len(self.dose))
            ]

        return headers, data

    def get_tbl_data_sds(self, disttype: constants.DistType):
        if disttype == constants.DistType.log_normal:
            headers = "Dose|N|Sample SD|Approximate Sample GSD|Model Fitted GSD|Scaled Residual"
            data = [
                (
                    self.dose[idx],
                    self.size[idx],
                    self.obs_sd[idx],
                    self.calc_sd[idx],
                    self.est_sd[idx],
                    self.residual[idx],
                )
                for idx in range(len(self.dose))
            ]
        else:
            headers = "Dose|N|Sample SD|Model Fitted SD"
            data = [
                (
                    self.dose[idx],
                    self.size[idx],
                    self.obs_sd[idx],
                    self.est_sd[idx],
                )
                for idx in range(len(self.dose))
            ]

        return headers, data

    def tbl(self, disttype: constants.DistType) -> str:
        mean_headers, mean_data = self.get_tbl_data_means(disttype)
        sd_headers, sd_data = self.get_tbl_data_sds(disttype)
        return "\n".join(
            [
                pretty_table(mean_data, mean_headers.split("|")),
                pretty_table(sd_data, sd_headers.split("|")),
            ]
        )

    def n(self) -> int:
        return self.dose.size


class ContinuousDeviance(BaseModel):
    names: list[str]
    loglikelihoods: list[float]
    num_params: list[int]
    aics: list[float]

    @classmethod
    def from_model(cls, model) -> Self:
        aod = model.structs.result.aod
        return cls(
            names=["A1", "A2", "A3", "fitted", "reduced"],
            loglikelihoods=aod.LL,
            num_params=aod.nParms,
            aics=aod.AIC,
        )

    def tbl(self) -> str:
        headers = "Model|Log-Likelihood|# Params|AIC".split("|")
        data = []
        for name, loglikelihood, num_param, aic in zip(
            self.names, self.loglikelihoods, self.num_params, self.aics, strict=True
        ):
            data.append([name, loglikelihood, num_param, aic])
        return pretty_table(data, headers)


class ContinuousTests(BaseModel):
    names: list[str]
    ll_ratios: list[float]
    dfs: list[float]
    p_values: list[float]

    @classmethod
    def from_model(cls, model) -> Self:
        tests = model.structs.result.aod.TOI
        return cls(
            names=["Test 1", "Test 2", "Test 3", "Test 4"],
            ll_ratios=tests.llRatio,
            dfs=tests.DF,
            p_values=tests.pVal,
        )

    def tbl(self) -> str:
        headers = "Name|-2 * Log(Likelihood Ratio)|Test d.f.|P-Value".split("|")
        data = []
        for name, ll_ratio, df, p_value in zip(
            self.names, self.ll_ratios, self.dfs, self.p_values, strict=True
        ):
            data.append([name, ll_ratio, df, p_value])
        return pretty_table(data, headers)


class ContinuousPlotting(BaseModel):
    dr_x: NumpyFloatArray
    dr_y: NumpyFloatArray
    bmdl_y: float
    bmd_y: float
    bmdu_y: float

    @classmethod
    def from_model(cls, model, params) -> Self:
        summary = model.structs.result.bmdsRes
        xs = np.array([summary.BMDL, summary.BMD, summary.BMDU])
        extra_values = [summary.BMD] if summary.BMD >= 0 else []
        dr_x = model.dataset.dose_linspace(extra_values=extra_values)
        dr_y = clean_array(model.dr_curve(dr_x, params))
        with warnings.catch_warnings():
            # xs may have a nan if a BMDU is not calculated; filter power warnings
            warnings.filterwarnings("ignore", message=".*invalid value encountered in power.*")
            critical_ys = clean_array(model.dr_curve(xs, params))
        critical_ys[critical_ys <= 0] = constants.BMDS_BLANK_VALUE
        return cls(
            dr_x=dr_x,
            dr_y=dr_y,
            bmdl_y=critical_ys[0],
            bmd_y=critical_ys[1],
            bmdu_y=critical_ys[2],
        )


class ContinuousResult(BaseModel):
    bmdl: float
    bmd: float
    bmdu: float
    has_completed: bool
    fit: ContinuousModelResult
    gof: ContinuousGof
    parameters: ContinuousParameters
    deviance: ContinuousDeviance
    tests: ContinuousTests
    plotting: ContinuousPlotting

    def tbl(self) -> str:
        data = [
            ["BMD", self.bmd],
            ["BMDL", self.bmdl],
            ["BMDU", self.bmdu],
            ["AIC", self.fit.aic],
            ["Log-Likelihood", self.fit.loglikelihood],
            ["P-Value", self.tests.p_values[3]],
            ["Model d.f.", self.tests.dfs[3]],
        ]
        return pretty_table(data, "")

    def text(self, dataset: ContinuousDatasets, settings: ContinuousModelSettings) -> str:
        return multi_lstrip(
            f"""
        Modeling Summary:
        {self.tbl()}

        Model Parameters:
        {self.parameters.tbl()}

        Goodness of Fit:
        {self.gof.tbl(disttype=settings.disttype)}

        Likelihoods:
        {self.deviance.tbl()}

        Tests of Mean and Variance Fits:
        {self.tests.tbl()}
        {CONTINUOUS_TEST_FOOTNOTES}
        """
        )

    @classmethod
    def from_model(cls, model) -> Self:
        summary = model.structs.result.bmdsRes
        params = ContinuousParameters.from_model(model)
        return cls(
            bmdl=summary.BMDL,
            bmd=summary.BMD,
            bmdu=summary.BMDU,
            has_completed=summary.validResult,
            fit=ContinuousModelResult.from_model(model),
            gof=ContinuousGof.from_model(model),
            parameters=params,
            deviance=ContinuousDeviance.from_model(model),
            tests=ContinuousTests.from_model(model),
            plotting=ContinuousPlotting.from_model(model, params.values),
        )

    def update_record(self, d: dict) -> None:
        """Update data record for a tabular-friendly export"""
        d.update(
            bmdl=self.bmdl,
            bmd=self.bmd,
            bmdu=self.bmdu,
            aic=self.fit.aic,
            loglikelihood=self.fit.loglikelihood,
            p_value1=self.tests.p_values[0],
            p_value2=self.tests.p_values[1],
            p_value3=self.tests.p_values[2],
            p_value4=self.tests.p_values[3],
            model_dof=self.tests.dfs[3],
            residual_of_interest=self.gof.roi,
            residual_at_lowest_dose=self.gof.residual[0],
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
                return self.fit.aic
            case "dof":
                return self.tests.dfs[3]
            case "pvalue":
                return self.tests.p_values[3]
            case "roi":
                return self.gof.roi
            case "roi_control":
                return self.gof.residual[0]
            case "n_params":
                return len(self.parameters.values)
            case _:  # pragma: no cover
                raise ValueError(f"Unknown parameter: {parameter}")
