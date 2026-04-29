import warnings
from typing import Self

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field

from .. import bmdscore, constants
from ..constants import BMDS_BLANK_VALUE, PriorClass
from ..models.dichotomous import BmdModelDichotomous
from .common import clean_array, inspect_cpp_obj
from .continuous import NumpyFloatArray
from .dichotomous import (
    DichotomousAnalysisOfDeviance,
    DichotomousModelResult,
    DichotomousParameters,
    DichotomousPgofResult,
    DichotomousPlotting,
    DichotomousResult,
)


class ModelAverageResult(BaseModel):
    pass


class ModelAveragePerModelSummary(BaseModel):
    bmdl: float
    bmd: float
    bmdu: float
    prior: float
    posterior: float


class DichotomousModelAverage:
    def __init__(
        self,
        dataset,
        models: list[BmdModelDichotomous],
        model_weights: npt.NDArray,
        prior_class: PriorClass = PriorClass.bayesian,
        weight_option: int | None = None,
    ):
        first = models[0].structs.analysis
        analysis = bmdscore.python_dichotomous_analysis()
        analysis.BMD_type = first.BMD_type
        analysis.BMR = first.BMR
        analysis.alpha = first.alpha
        analysis.samples = first.samples
        analysis.burnin = first.burnin
        analysis.degree = first.degree
        analysis.Y = dataset.incidences
        analysis.n_group = dataset.ns
        analysis.doses = dataset.doses
        analysis.n = dataset.num_dose_groups

        average = bmdscore.python_dichotomousMA_analysis()
        average.nmodels = len(models)
        average.nparms = [model.structs.result.nparms for model in models]
        average.actual_parms = [model.structs.result.nparms for model in models]
        average.prior_cols = [model.structs.analysis.prior_cols for model in models]
        average.models = [model.structs.analysis.model for model in models]
        average.priors = [model.structs.analysis.prior for model in models]
        average.modelPriors = model_weights
        average.pyDA = analysis
        average.weightOption = 1 if weight_option is None else int(weight_option)
        average.datatype = (
            bmdscore.loud_datatype.l_dichotomous.value
            if prior_class is PriorClass.bayesian_loud
            else 0
        )

        bmdsRes = bmdscore.BMDSMA_results()
        bmdsRes.BMD = np.full(average.nmodels, -9999)
        bmdsRes.BMDL = np.full(average.nmodels, -9999)
        bmdsRes.BMDU = np.full(average.nmodels, -9999)
        bmdsRes.ebUpper = np.full(analysis.n, -9999)
        bmdsRes.ebLower = np.full(analysis.n, -9999)

        result = bmdscore.python_dichotomousMA_result()
        result.nmodels = len(models)
        result.dist_numE = 200
        result.post_probs = [0.0] * result.nmodels
        if prior_class is PriorClass.bayesian_loud:
            result.bmd_dist = [0.0] * int(analysis.samples)
            result.models = [
                bmdscore.python_dichotomous_model_result() for _ in range(result.nmodels)
            ]
            for i, model in enumerate(models):
                nparms = int(model.structs.result.nparms)
                iters = int(analysis.samples)
                loud = result.models[i].loudRes
                loud.parms = [0.0] * (nparms * iters)
                loud.BMD = [0.0] * iters
        else:
            result.bmd_dist = [0.0] * (2 * result.dist_numE)
            result.models = [model.structs.result for model in models]
        result.bmdsRes = bmdsRes

        self.analysis = analysis
        self.average = average
        self.result = result
        self.bmdsRes = result.bmdsRes  # use this version; copied on assignment above

    def execute(self) -> "DichotomousModelAverageResult":
        if self.average.datatype == bmdscore.loud_datatype.l_dichotomous.value:
            bmdscore.pythonBMDSLoud(self.average, self.result)
        else:
            bmdscore.pythonBMDSDichoMA(self.average, self.result)

    def __str__(self) -> str:
        return "\n".join(
            [
                inspect_cpp_obj(self.analysis),
                inspect_cpp_obj(self.result),
            ]
        )


class DichotomousModelAverageResult(ModelAverageResult):
    """
    Model average fit
    """

    bmdl: float
    bmd: float
    bmdu: float
    bmdl_y: float
    bmd_y: float
    bmdu_y: float
    bmd_dist: NumpyFloatArray
    priors: NumpyFloatArray
    posteriors: NumpyFloatArray
    model_bmd_dist: list[NumpyFloatArray] = Field(default_factory=list)
    model_parm_dist: list[NumpyFloatArray] = Field(default_factory=list)
    dr_x: NumpyFloatArray
    dr_y: NumpyFloatArray

    @staticmethod
    def _param_names(model, n_params: int) -> list[str]:
        names = model.get_param_names()
        if len(names) >= n_params:
            return names[:n_params]
        return [*names, *[f"param_{i + 1}" for i in range(len(names), n_params)]]

    @staticmethod
    def _safe_dr_curve(model, dr_x, params) -> np.ndarray | None:
        try:
            return clean_array(model.dr_curve(dr_x, params))
        except (FloatingPointError, IndexError, ValueError):
            return None

    @staticmethod
    def _resize_draws(draws: np.ndarray, n_draws: int) -> np.ndarray:
        draws = np.asarray(draws, dtype=float)
        if draws.ndim == 1:
            out = np.full(n_draws, np.nan, dtype=float)
            out[: min(n_draws, draws.size)] = draws[:n_draws]
            return out

        out = np.full((n_draws, draws.shape[1]), np.nan, dtype=float)
        rows = min(n_draws, draws.shape[0])
        out[:rows, :] = draws[:rows, :]
        return out

    @classmethod
    def from_cpp(cls, analysis: DichotomousModelAverage, model_results) -> Self:
        # only keep positive finite values
        is_loud = analysis.average.datatype == bmdscore.loud_datatype.l_dichotomous.value
        if is_loud:
            arr = np.asarray(analysis.result.bmd_dist, dtype=float)
            arr = arr[np.isfinite(arr)]
        else:
            arr = np.array(analysis.result.bmd_dist).reshape(2, analysis.result.dist_numE).T
            arr = arr[np.isfinite(arr[:, 0])]
            arr = arr[arr[:, 0] > 0]

        # calculate dr_y for model averaging
        priors = np.array(analysis.average.modelPriors)
        posteriors = np.array(analysis.result.post_probs)
        model_bmd: list[np.ndarray] = []
        model_parms: list[np.ndarray] = []
        if is_loud:
            n_draws = arr.size
            for result in analysis.result.models:
                bmd_draws = np.asarray(result.loudRes.BMD, dtype=float)
                n_draw = bmd_draws.size
                model_bmd.append(cls._resize_draws(bmd_draws, n_draws))

                parm_draws = np.asarray(result.loudRes.parms, dtype=float)
                if parm_draws.ndim == 1 and n_draw > 0 and parm_draws.size % n_draw == 0:
                    parm_draws = parm_draws.reshape(n_draw, parm_draws.size // n_draw)
                elif (
                    parm_draws.ndim == 2
                    and n_draw > 0
                    and parm_draws.shape[0] != n_draw
                    and parm_draws.shape[1] == n_draw
                ):
                    parm_draws = parm_draws.T
                model_parms.append(cls._resize_draws(parm_draws, n_draws))

        if is_loud:
            models = model_results
            dr_x = np.asarray(models[0].dataset.dose_linspace(), dtype=float)
            values = []
            valid_posteriors = []
            for idx, (model, parm_draws) in enumerate(zip(models, model_parms, strict=True)):
                finite_parm_draws = parm_draws[np.isfinite(parm_draws).all(axis=1)]
                params = (
                    np.nanmedian(finite_parm_draws, axis=0)
                    if finite_parm_draws.size
                    else np.array([])
                )
                dr_y_ = cls._safe_dr_curve(model, dr_x, params)
                if dr_y_ is not None:
                    values.append(dr_y_)
                    valid_posteriors.append(posteriors[idx])
            values = np.asarray(values, dtype=float)
            valid_posteriors = np.asarray(valid_posteriors, dtype=float)
            if values.size and valid_posteriors.sum() > 0:
                dr_y = (valid_posteriors / valid_posteriors.sum()) @ values
            else:
                dr_y = np.full_like(dr_x, BMDS_BLANK_VALUE, dtype=float)
        else:
            values = np.array([result.plotting.dr_y for result in model_results])
            dr_x = model_results[0].plotting.dr_x
            dr_y = values.T.dot(posteriors)
        bmds = [analysis.bmdsRes.BMDL_MA, analysis.bmdsRes.BMD_MA, analysis.bmdsRes.BMDU_MA]
        bmds_ys = np.interp(bmds, dr_x, dr_y)
        return cls(
            bmdl=bmds[0],
            bmd=bmds[1],
            bmdu=bmds[2],
            bmdl_y=bmds_ys[0],
            bmd_y=bmds_ys[1],
            bmdu_y=bmds_ys[2],
            bmd_dist=arr if is_loud else arr.T,
            priors=priors,
            posteriors=posteriors,
            model_bmd_dist=model_bmd,
            model_parm_dist=model_parms,
            dr_x=dr_x,
            dr_y=dr_y,
        )

    def update_record(self, d: dict) -> None:
        """Update data record for a tabular-friendly export"""
        d.update(
            bmdl=self.bmdl,
            bmd=self.bmd,
            bmdu=self.bmdu,
        )

    def update_record_weights(self, d: dict, index: int) -> None:
        """Update data record for a tabular-friendly export"""
        d.update(
            model_prior=self.priors[index],
            model_posterior=self.posteriors[index],
        )

    def model_summary(self, index: int, alpha: float) -> ModelAveragePerModelSummary:
        draws = np.asarray(self.model_bmd_dist[index], dtype=float)
        draws = draws[np.isfinite(draws)]
        prior = float(self.priors[index])
        posterior = float(self.posteriors[index])
        if draws.size == 0:
            return ModelAveragePerModelSummary(
                bmdl=BMDS_BLANK_VALUE,
                bmd=BMDS_BLANK_VALUE,
                bmdu=BMDS_BLANK_VALUE,
                prior=prior,
                posterior=posterior,
            )
        return ModelAveragePerModelSummary(
            bmdl=float(np.quantile(draws, alpha)),
            bmd=float(np.quantile(draws, 0.5)),
            bmdu=float(np.quantile(draws, 1 - alpha)),
            prior=prior,
            posterior=posterior,
        )

    def sync_model_result(self, model, index: int) -> None:
        summary = self.model_summary(index, model.settings.alpha)
        draws = np.asarray(self.model_bmd_dist[index], dtype=float)
        draws = draws[np.isfinite(draws)]
        raw_parm_draws = np.asarray(self.model_parm_dist[index], dtype=float)
        if raw_parm_draws.ndim == 2:
            n_params = raw_parm_draws.shape[1]
            parm_draws = raw_parm_draws[np.isfinite(raw_parm_draws).all(axis=1)]
        else:
            n_params = raw_parm_draws.size
            parm_draws = raw_parm_draws[np.isfinite(raw_parm_draws)]

        if parm_draws.ndim == 2 and parm_draws.size:
            param_values = np.nanmedian(parm_draws, axis=0)
            lower_ci = np.nanquantile(parm_draws, model.settings.alpha, axis=0)
            upper_ci = np.nanquantile(parm_draws, 1 - model.settings.alpha, axis=0)
            se = np.nanstd(parm_draws, axis=0)
            cov = np.cov(parm_draws, rowvar=False)
        else:
            param_values = np.full(n_params, BMDS_BLANK_VALUE, dtype=float)
            lower_ci = np.full_like(param_values, BMDS_BLANK_VALUE, dtype=float)
            upper_ci = np.full_like(param_values, BMDS_BLANK_VALUE, dtype=float)
            se = np.full_like(param_values, BMDS_BLANK_VALUE, dtype=float)
            cov = np.full((param_values.size, param_values.size), BMDS_BLANK_VALUE, dtype=float)

        priors = DichotomousParameters.get_priors(model)
        param_names = self._param_names(model, param_values.size)
        parameters = DichotomousParameters(
            names=param_names,
            values=param_values,
            se=se,
            lower_ci=lower_ci,
            upper_ci=upper_ci,
            bounded=np.zeros(param_values.size, dtype=float),
            cov=np.asarray(cov, dtype=float).reshape(param_values.size, param_values.size),
            prior_type=priors[0][: param_values.size],
            prior_initial_value=priors[1][: param_values.size],
            prior_stdev=priors[2][: param_values.size],
            prior_min_value=priors[3][: param_values.size],
            prior_max_value=priors[4][: param_values.size],
        )

        fit = DichotomousModelResult(
            loglikelihood=BMDS_BLANK_VALUE,
            aic=BMDS_BLANK_VALUE,
            bic_equiv=BMDS_BLANK_VALUE,
            chisq=BMDS_BLANK_VALUE,
            model_df=BMDS_BLANK_VALUE,
            total_df=BMDS_BLANK_VALUE,
            bmd_dist=draws,
        )

        n = model.dataset.num_dose_groups
        gof = DichotomousPgofResult(
            expected=[BMDS_BLANK_VALUE] * n,
            residual=[BMDS_BLANK_VALUE] * n,
            eb_lower=[BMDS_BLANK_VALUE] * n,
            eb_upper=[BMDS_BLANK_VALUE] * n,
            test_statistic=BMDS_BLANK_VALUE,
            p_value=getattr(model.structs.result.loudRes, "pval", BMDS_BLANK_VALUE),
            roi=BMDS_BLANK_VALUE,
            df=BMDS_BLANK_VALUE,
        )

        deviance = DichotomousAnalysisOfDeviance(
            names=["Full model", "Fitted model", "Reduced model"],
            ll=[BMDS_BLANK_VALUE] * 3,
            params=[BMDS_BLANK_VALUE] * 3,
            deviance=[BMDS_BLANK_VALUE] * 3,
            df=[BMDS_BLANK_VALUE] * 3,
            p_value=[BMDS_BLANK_VALUE] * 3,
        )

        extra_values = [summary.bmd] if summary.bmd >= 0 else []
        dr_x = model.dataset.dose_linspace(extra_values=extra_values)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*invalid value encountered.*")
            warnings.filterwarnings("ignore", message=".*divide by zero encountered.*")
            dr_y = self._safe_dr_curve(model, dr_x, parameters.values)
        if dr_y is None:
            dr_y = np.full_like(np.asarray(dr_x, dtype=float), BMDS_BLANK_VALUE, dtype=float)

        xs = np.asarray([summary.bmdl, summary.bmd, summary.bmdu], dtype=float)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*invalid value encountered.*")
            warnings.filterwarnings("ignore", message=".*divide by zero encountered.*")
            critical_ys = self._safe_dr_curve(model, xs, parameters.values)
        if critical_ys is None:
            critical_ys = np.full_like(xs, BMDS_BLANK_VALUE, dtype=float)
        critical_ys[critical_ys <= 0] = constants.BMDS_BLANK_VALUE

        plotting = DichotomousPlotting(
            dr_x=dr_x,
            dr_y=dr_y,
            bmdl_y=float(critical_ys[0]),
            bmd_y=float(critical_ys[1]),
            bmdu_y=float(critical_ys[2]),
        )

        model.results = DichotomousResult(
            has_completed=True,
            bmdl=summary.bmdl,
            bmd=summary.bmd,
            bmdu=summary.bmdu,
            slope_factor=None,
            fit=fit,
            gof=gof,
            parameters=parameters,
            deviance=deviance,
            plotting=plotting,
        )
