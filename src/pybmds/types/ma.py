import warnings
from typing import Self

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, model_serializer

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
        n_chains = models[0].settings.n_chains
        seed = models[0].settings.seed
        if hasattr(analysis, "chains"):
            analysis.chains = n_chains
        if seed is not None and hasattr(analysis, "seed"):
            analysis.seed = seed

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
        if seed is not None and hasattr(average, "seed"):
            average.seed = seed

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
            result.models = [
                bmdscore.python_dichotomous_model_result() for _ in range(result.nmodels)
            ]
        else:
            result.bmd_dist = [0.0] * (2 * result.dist_numE)
            result.models = [model.structs.result for model in models]
        result.bmdsRes = bmdsRes

        self.analysis = analysis
        self.average = average
        self.result = result
        self.n_chains = n_chains
        self.seed = seed
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
    model_p_values: list[float] = Field(default_factory=list)
    model_waics: list[float] = Field(default_factory=list)
    model_int_factors: list[float] = Field(default_factory=list)
    model_loglikelihoods: list[float] = Field(default_factory=list)
    dr_x: NumpyFloatArray
    dr_y: NumpyFloatArray

    def without_loud_draws(self) -> Self:
        return self.model_copy(
            update={
                "bmd_dist": np.empty((0,), dtype=float),
                "model_bmd_dist": [],
                "model_parm_dist": [],
            }
        )

    @staticmethod
    def _json_safe_draws(draws) -> list:
        arr = np.asarray(draws, dtype=float)
        values = arr.astype(object)
        values[~np.isfinite(arr)] = None
        return values.tolist()

    @staticmethod
    def _apply_paired_valid_draw_mask(
        bmd_draws: np.ndarray, parm_draws: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray]:
        bmd = np.asarray(bmd_draws, dtype=float).copy()
        parms = np.asarray(parm_draws, dtype=float).copy()

        if bmd.size == 0 or parms.size == 0:
            return bmd, parms

        if parms.ndim == bmd.ndim + 1 and parms.shape[:-1] == bmd.shape:
            valid = np.isfinite(bmd) & np.isfinite(parms).all(axis=-1)
        elif bmd.ndim == 1 and parms.ndim == 2 and parms.shape[0] == bmd.shape[0]:
            valid = np.isfinite(bmd) & np.isfinite(parms).all(axis=1)
        elif parms.ndim == 2 and parms.shape[0] == bmd.size:
            valid = (np.isfinite(bmd.reshape(-1)) & np.isfinite(parms).all(axis=1)).reshape(
                bmd.shape
            )
        else:
            return bmd, parms

        bmd[~valid] = np.nan
        if parms.ndim == 2 and parms.shape[0] == bmd.size and valid.shape != parms.shape[:1]:
            parms[~valid.reshape(-1), :] = np.nan
        else:
            parms[~valid, ...] = np.nan
        return bmd, parms

    @model_serializer(mode="wrap")
    def serialize_model(self, handler):
        data = handler(self)

        bmd_dist = np.asarray(self.bmd_dist, dtype=float)
        if bmd_dist.ndim == 1:
            data["bmd_dist"] = self._json_safe_draws(self.bmd_dist)

        if self.model_bmd_dist and self.model_parm_dist:
            data["model_bmd_dist"] = [self._json_safe_draws(draws) for draws in self.model_bmd_dist]
            data["model_parm_dist"] = [
                self._json_safe_draws(draws) for draws in self.model_parm_dist
            ]

        return data

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

    @staticmethod
    def _loud_weights(scores: list[float]) -> np.ndarray:
        scores = np.asarray(scores, dtype=float)
        valid = np.isfinite(scores) & (scores != BMDS_BLANK_VALUE)
        weights = np.zeros(scores.size, dtype=float)
        if not valid.any():
            return np.full(scores.size, BMDS_BLANK_VALUE, dtype=float)
        finite_scores = scores[valid]
        weights[valid] = np.exp(finite_scores - finite_scores.max())
        weights_sum = weights.sum()
        if weights_sum <= 0 or not np.isfinite(weights_sum):
            return np.full(scores.size, BMDS_BLANK_VALUE, dtype=float)
        return weights / weights_sum

    @staticmethod
    def _reshape_flat_draws(draws: np.ndarray, n_chains: int) -> np.ndarray:
        draws = np.asarray(draws, dtype=float)
        if draws.ndim != 1 or n_chains <= 1:
            return draws
        n = draws.shape[0]
        if n % n_chains != 0:
            return draws
        return draws.reshape(n_chains, n // n_chains)

    @staticmethod
    def _combined_loud_result(result):
        combined = getattr(result, "combinedLoudRes", None)
        if combined is not None and np.asarray(getattr(combined, "BMD", []), dtype=float).size:
            return combined

        loud_res = getattr(result, "loudRes", None)
        if isinstance(loud_res, list) and loud_res:
            return loud_res
        return loud_res

    @classmethod
    def _loud_draws(cls, loud, n_chains: int) -> tuple[np.ndarray, np.ndarray]:
        if isinstance(loud, list):
            bmd_by_chain = [np.asarray(chain.BMD, dtype=float) for chain in loud]
            parm_by_chain = []
            for chain, bmd in zip(loud, bmd_by_chain, strict=True):
                parms = np.asarray(chain.parms, dtype=float)
                if parms.ndim == 1 and bmd.size > 0 and parms.size % bmd.size == 0:
                    parms = parms.reshape(bmd.size, parms.size // bmd.size)
                elif (
                    parms.ndim == 2
                    and bmd.size > 0
                    and parms.shape[0] != bmd.size
                    and parms.shape[1] == bmd.size
                ):
                    parms = parms.T
                parm_by_chain.append(parms)
            return np.asarray(bmd_by_chain, dtype=float), np.asarray(parm_by_chain, dtype=float)

        bmd_draws = np.asarray(getattr(loud, "BMD", []), dtype=float)
        n_draw = bmd_draws.size
        parm_draws = np.asarray(getattr(loud, "parms", []), dtype=float)
        if parm_draws.ndim == 1 and n_draw > 0 and parm_draws.size % n_draw == 0:
            parm_draws = parm_draws.reshape(n_draw, parm_draws.size // n_draw)
        elif (
            parm_draws.ndim == 2
            and n_draw > 0
            and parm_draws.shape[0] != n_draw
            and parm_draws.shape[1] == n_draw
        ):
            parm_draws = parm_draws.T
        return cls._reshape_flat_draws(bmd_draws, n_chains), cls._reshape_flat_draws(
            parm_draws, n_chains
        )

    @classmethod
    def _loud_posteriors(
        cls,
        posteriors: np.ndarray,
        weight_option: int,
        waics: list[float],
        int_factors: list[float],
    ) -> np.ndarray:
        valid = (
            np.isfinite(posteriors)
            & (posteriors != BMDS_BLANK_VALUE)
            & (posteriors >= 0)
            & (posteriors <= 1)
        )
        if valid.any() and np.isclose(posteriors[valid].sum(), 1.0):
            return posteriors

        if weight_option == 1:
            return cls._loud_weights(waics)
        if weight_option == 2:
            return cls._loud_weights(int_factors)

        waic_weights = cls._loud_weights(waics)
        int_factor_weights = cls._loud_weights(int_factors)
        usable = [
            weights
            for weights in (waic_weights, int_factor_weights)
            if np.isfinite(weights).any() and (weights != BMDS_BLANK_VALUE).any()
        ]
        if not usable:
            return np.full(posteriors.size, BMDS_BLANK_VALUE, dtype=float)
        posteriors = np.mean(usable, axis=0)
        return posteriors / posteriors.sum()

    @classmethod
    def from_cpp(cls, analysis: DichotomousModelAverage, model_results) -> Self:
        # only keep positive finite values
        is_loud = analysis.average.datatype == bmdscore.loud_datatype.l_dichotomous.value
        if is_loud:
            n_chains = int(getattr(analysis, "n_chains", 1) or 1)
            arr = cls._reshape_flat_draws(
                np.asarray(analysis.result.bmd_dist, dtype=float), n_chains
            )
        else:
            arr = np.array(analysis.result.bmd_dist).reshape(2, analysis.result.dist_numE).T
            arr = arr[np.isfinite(arr[:, 0])]
            arr = arr[arr[:, 0] > 0]

        # calculate dr_y for model averaging
        priors = np.array(analysis.average.modelPriors)
        posteriors = np.array(analysis.result.post_probs)
        model_bmd: list[np.ndarray] = []
        model_parms: list[np.ndarray] = []
        model_p_values: list[float] = []
        model_waics: list[float] = []
        model_int_factors: list[float] = []
        model_loglikelihoods: list[float] = []
        if is_loud:
            for result in analysis.result.models:
                loud = cls._combined_loud_result(result)
                model_p_values.append(float(getattr(loud, "pval", BMDS_BLANK_VALUE)))
                model_waics.append(float(getattr(loud, "waic", BMDS_BLANK_VALUE)))
                model_int_factors.append(float(getattr(loud, "int_factor", BMDS_BLANK_VALUE)))
                model_loglikelihoods.append(float(getattr(loud, "ll", BMDS_BLANK_VALUE)))
                bmd_draws, parm_draws = cls._loud_draws(loud, n_chains)
                bmd_draws, parm_draws = cls._apply_paired_valid_draw_mask(bmd_draws, parm_draws)
                model_bmd.append(bmd_draws)
                model_parms.append(parm_draws)

        if is_loud:
            posteriors = cls._loud_posteriors(
                posteriors, analysis.average.weightOption, model_waics, model_int_factors
            )
            models = model_results
            dr_x = np.asarray(models[0].dataset.dose_linspace(), dtype=float)
            values = []
            valid_posteriors = []
            for idx, (model, parm_draws) in enumerate(zip(models, model_parms, strict=True)):
                finite_parm_draws = parm_draws.reshape(-1, parm_draws.shape[-1])
                finite_parm_draws = finite_parm_draws[np.isfinite(finite_parm_draws).all(axis=1)]
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
            model_p_values=model_p_values,
            model_waics=model_waics,
            model_int_factors=model_int_factors,
            model_loglikelihoods=model_loglikelihoods,
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

    def sync_model_result(self, model, index: int, cpp_result=None) -> None:
        summary = self.model_summary(index, model.settings.alpha)
        draws = np.asarray(self.model_bmd_dist[index], dtype=float)
        draws = draws[np.isfinite(draws)]
        raw_parm_draws = np.asarray(self.model_parm_dist[index], dtype=float)
        if raw_parm_draws.ndim >= 2:
            n_params = raw_parm_draws.shape[-1]
            flat_parm_draws = raw_parm_draws.reshape(-1, n_params)
            parm_draws = flat_parm_draws[np.isfinite(flat_parm_draws).all(axis=1)]
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
            loglikelihood=(
                self.model_loglikelihoods[index]
                if index < len(self.model_loglikelihoods)
                else BMDS_BLANK_VALUE
            ),
            aic=BMDS_BLANK_VALUE,
            bic_equiv=BMDS_BLANK_VALUE,
            chisq=BMDS_BLANK_VALUE,
            model_df=BMDS_BLANK_VALUE,
            total_df=BMDS_BLANK_VALUE,
            bmd_dist=draws,
        )

        if cpp_result is not None:
            gof = DichotomousPgofResult.from_cpp(cpp_result.gof, summary.bmd, model.dataset.doses)
            gof.p_value = (
                self.model_p_values[index] if index < len(self.model_p_values) else BMDS_BLANK_VALUE
            )
        else:
            n = model.dataset.num_dose_groups
            gof = DichotomousPgofResult(
                expected=[BMDS_BLANK_VALUE] * n,
                residual=[BMDS_BLANK_VALUE] * n,
                eb_lower=[BMDS_BLANK_VALUE] * n,
                eb_upper=[BMDS_BLANK_VALUE] * n,
                test_statistic=BMDS_BLANK_VALUE,
                p_value=(
                    self.model_p_values[index]
                    if index < len(self.model_p_values)
                    else BMDS_BLANK_VALUE
                ),
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
            summary_p_value=(
                self.model_p_values[index] if index < len(self.model_p_values) else BMDS_BLANK_VALUE
            ),
            summary_waic=(
                self.model_waics[index] if index < len(self.model_waics) else BMDS_BLANK_VALUE
            ),
        )
