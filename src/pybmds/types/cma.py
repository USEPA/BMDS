import warnings
from typing import Self

import numpy as np
import numpy.typing as npt
from pydantic import BaseModel, Field, model_serializer

from .. import bmdscore, constants
from ..models.continuous import BmdModelContinuous
from .common import clean_array, inspect_cpp_obj
from .continuous import (
    ContinuousDeviance,
    ContinuousGof,
    ContinuousModelResult,
    ContinuousParameters,
    ContinuousPlotting,
    ContinuousResult,
    ContinuousTests,
    NumpyFloatArray,
)


class ModelAverageResult(BaseModel):
    pass


class ModelAveragePerModelSummary(BaseModel):
    bmdl: float
    bmd: float
    bmdu: float
    prior: float
    posterior: float


class ContinuousModelAverage:
    def __init__(
        self,
        dataset,
        models: list[BmdModelContinuous],
        model_weights: npt.NDArray,
        weight_option: int | None = None,
    ):
        first = models[0].structs.analysis
        analysis = bmdscore.python_continuous_analysis()
        analysis.BMD_type = first.BMD_type
        analysis.BMR = first.BMR
        analysis.alpha = first.alpha
        analysis.suff_stat = first.suff_stat
        analysis.disttype = first.disttype
        analysis.isIncreasing = first.isIncreasing
        analysis.detectAdvDir = first.detectAdvDir
        analysis.samples = first.samples
        analysis.burnin = first.burnin
        if analysis.suff_stat is True:
            analysis.Y = dataset.means
            analysis.sd = dataset.stdevs
            analysis.n_group = dataset.ns
            analysis.doses = dataset.doses
            analysis.n = dataset.num_dose_groups
        elif analysis.suff_stat is False:
            analysis.Y = dataset.responses
            analysis.sd = []
            analysis.n_group = []
            analysis.doses = dataset.individual_doses
            analysis.n = len(dataset.individual_doses)
        else:
            raise ValueError(f"Unsupported dataset dtype: {dataset.dtype}")
        n_chains = models[0].settings.n_chains
        seed = models[0].settings.seed
        if hasattr(analysis, "chains"):
            analysis.chains = n_chains
        if seed is not None and hasattr(analysis, "seed"):
            analysis.seed = seed

        average = bmdscore.python_continuousMA_analysis()
        average.nmodels = len(models)
        average.nparms = [model.structs.result.nparms for model in models]
        average.actual_parms = [model.structs.result.nparms for model in models]
        average.prior_cols = [model.structs.analysis.prior_cols for model in models]
        average.models = [int(m.structs.analysis.model) for m in models]
        average.disttype = [int(m.structs.analysis.disttype) for m in models]
        average.weightOption = 1 if weight_option is None else int(weight_option)

        if analysis.suff_stat is True:
            average.datatype = bmdscore.loud_datatype.l_summary
        elif analysis.suff_stat is False:
            average.datatype = bmdscore.loud_datatype.l_individual

        average.priors = [model.structs.analysis.prior for model in models]

        average.modelPriors = [float(x) for x in model_weights]
        average.pyCA = analysis
        if seed is not None and hasattr(average, "seed"):
            average.seed = seed

        bmdsRes = bmdscore.BMDSMA_results()
        bmdsRes.BMD_MA = -9999.0
        bmdsRes.BMDL_MA = -9999.0
        bmdsRes.BMDU_MA = -9999.0
        bmdsRes.BMD = np.full(len(models), -9999.0)
        bmdsRes.BMDL = np.full(len(models), -9999.0)
        bmdsRes.BMDU = np.full(len(models), -9999.0)
        bmdsRes.ebUpper = np.full(analysis.n, -9999)
        bmdsRes.ebLower = np.full(analysis.n, -9999)

        result = bmdscore.python_continuousMA_result()
        result.nmodels = len(models)
        result.dist_numE = 200

        result.post_probs = [0.0] * result.nmodels
        result.models = [bmdscore.python_continuous_model_result() for _ in range(result.nmodels)]

        result.bmdsRes = bmdsRes

        self.analysis = analysis
        self.average = average
        self.result = result
        self.n_chains = n_chains
        self.seed = seed
        self.bmdsRes = result.bmdsRes  # use this version; copied on assignment above

    def execute(self) -> "ContinuousModelAverageResult":
        bmdscore.pythonBMDSLoud(self.average, self.result)
        return self

    def __str__(self) -> str:
        return "\n".join(
            [
                inspect_cpp_obj(self.analysis),
                inspect_cpp_obj(self.result),
            ]
        )


class ContinuousModelAverageResult(ModelAverageResult):
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
    model_bmd_dist: list[NumpyFloatArray]
    model_parm_dist: list[NumpyFloatArray]
    model_p_values: list[float] = Field(default_factory=list)
    model_waics: list[float] = Field(default_factory=list)
    model_loglikelihoods: list[float] = Field(default_factory=list)
    dr_x: NumpyFloatArray
    dr_y: NumpyFloatArray

    @staticmethod
    def _json_safe_draws(draws) -> list:
        arr = np.asarray(draws, dtype=float)
        values = arr.astype(object)
        values[~np.isfinite(arr)] = None
        return values.tolist()

    @model_serializer(mode="wrap")
    def serialize_model(self, handler):
        data = handler(self)

        data["bmd_dist"] = self._json_safe_draws(self.bmd_dist)
        data["model_bmd_dist"] = [self._json_safe_draws(draws) for draws in self.model_bmd_dist]
        data["model_parm_dist"] = [self._json_safe_draws(draws) for draws in self.model_parm_dist]
        return data

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
    def _loud_weights(scores: list[float]) -> np.ndarray:
        scores = np.asarray(scores, dtype=float)
        valid = np.isfinite(scores) & (scores != constants.BMDS_BLANK_VALUE)
        weights = np.zeros(scores.size, dtype=float)
        if not valid.any():
            return np.full(scores.size, constants.BMDS_BLANK_VALUE, dtype=float)
        finite_scores = scores[valid]
        weights[valid] = np.exp(finite_scores - finite_scores.max())
        weights_sum = weights.sum()
        if weights_sum <= 0 or not np.isfinite(weights_sum):
            return np.full(scores.size, constants.BMDS_BLANK_VALUE, dtype=float)
        return weights / weights_sum

    @classmethod
    def _loud_posteriors(cls, posteriors: np.ndarray, waics: list[float]) -> np.ndarray:
        posteriors = np.asarray(posteriors, dtype=float)
        valid = (
            np.isfinite(posteriors)
            & (posteriors != constants.BMDS_BLANK_VALUE)
            & (posteriors >= 0)
            & (posteriors <= 1)
        )
        if valid.any() and np.isclose(posteriors[valid].sum(), 1.0):
            return posteriors
        return cls._loud_weights(waics)

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
    def from_cpp(cls, analysis: ContinuousModelAverage, model_results=None) -> Self:
        n_chains = int(getattr(analysis, "n_chains", 1) or 1)
        ma_bmd = cls._reshape_flat_draws(
            np.asarray(analysis.result.bmd_dist, dtype=float), n_chains
        )

        model_bmd: list[np.ndarray] = []
        model_parms: list[np.ndarray] = []
        model_p_values: list[float] = []
        model_waics: list[float] = []
        model_loglikelihoods: list[float] = []

        for m in analysis.result.models:
            lr = cls._combined_loud_result(m)
            model_p_values.append(float(getattr(lr, "pval", constants.BMDS_BLANK_VALUE)))
            model_waics.append(float(getattr(lr, "waic", constants.BMDS_BLANK_VALUE)))
            model_loglikelihoods.append(float(getattr(lr, "ll", constants.BMDS_BLANK_VALUE)))
            bmd_draws, parm_draws = cls._loud_draws(lr, n_chains)
            model_bmd.append(bmd_draws)
            model_parms.append(parm_draws)

        priors = np.asarray(analysis.average.modelPriors, dtype=float)
        posteriors = np.asarray(analysis.result.post_probs, dtype=float)
        if posteriors.size == 0:
            raise RuntimeError("C++ MA did not populate post_probs (empty).")
        posteriors = cls._loud_posteriors(posteriors, model_waics)

        bmdsRes = analysis.result.bmdsRes
        bmds = np.asarray([bmdsRes.BMDL_MA, bmdsRes.BMD_MA, bmdsRes.BMDU_MA], dtype=float)
        if model_results is not None:
            values = np.asarray([r.plotting.dr_y for r in model_results], dtype=float)
            dr_x = np.asarray(model_results[0].plotting.dr_x, dtype=float)
            if posteriors.size != values.shape[0]:
                raise RuntimeError(
                    f"post_probs ({posteriors.size}) != nmodels ({values.shape[0]})."
                )
            dr_y = posteriors @ values
        else:
            dr_x = np.asarray(analysis.analysis.doses, dtype=float)
            dr_y = np.full(dr_x.shape, constants.BMDS_BLANK_VALUE, dtype=float)
        bmds_ys = np.interp(bmds, dr_x, dr_y)

        return cls(
            bmdl=float(bmds[0]),
            bmd=float(bmds[1]),
            bmdu=float(bmds[2]),
            bmdl_y=float(bmds_ys[0]),
            bmd_y=float(bmds_ys[1]),
            bmdu_y=float(bmds_ys[2]),
            bmd_dist=ma_bmd,
            model_bmd_dist=model_bmd,
            model_parm_dist=model_parms,
            model_p_values=model_p_values,
            model_waics=model_waics,
            model_loglikelihoods=model_loglikelihoods,
            priors=priors,
            posteriors=posteriors,
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
                bmdl=-9999.0,
                bmd=-9999.0,
                bmdu=-9999.0,
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
        """
        Update a bayesian LOUD model's public result object so downstream callers
        can use the per-model outputs stored in the LOUD MA result container.
        """
        summary = self.model_summary(index, model.settings.alpha)
        draws = np.asarray(self.model_bmd_dist[index], dtype=float)
        draws = draws[np.isfinite(draws)]
        parm_draws = np.asarray(self.model_parm_dist[index], dtype=float)
        if parm_draws.size:
            parameters = ContinuousParameters.from_loud_draws(
                model, parm_draws.reshape(-1, parm_draws.shape[-1])
            )
        elif model.results is not None:
            parameters = model.results.parameters
        else:
            raise RuntimeError(f"LOUD did not return parameter draws for {model.name()}.")

        if model.results is None:
            n = len(model.dataset.doses)
            empty = np.full(n, constants.BMDS_BLANK_VALUE, dtype=float)
            fit = ContinuousModelResult(
                dist=int(model.settings.disttype.value),
                loglikelihood=constants.BMDS_BLANK_VALUE,
                aic=constants.BMDS_BLANK_VALUE,
                bic_equiv=constants.BMDS_BLANK_VALUE,
                chisq=constants.BMDS_BLANK_VALUE,
                model_df=constants.BMDS_BLANK_VALUE,
                total_df=constants.BMDS_BLANK_VALUE,
                bmd_dist=draws,
            )
            if cpp_result is not None:
                gof = ContinuousGof.from_cpp(cpp_result.gof, summary.bmd, model.dataset.doses)
            else:
                gof = ContinuousGof(
                    dose=np.asarray(model.dataset.doses, dtype=float),
                    size=np.asarray(model.dataset.ns, dtype=float),
                    est_mean=empty,
                    calc_mean=empty,
                    obs_mean=np.asarray(model.dataset.means, dtype=float),
                    est_sd=empty,
                    calc_sd=empty,
                    obs_sd=np.asarray(model.dataset.stdevs, dtype=float),
                    residual=empty,
                    eb_lower=empty,
                    eb_upper=empty,
                    roi=constants.BMDS_BLANK_VALUE,
                )
            deviance = ContinuousDeviance(names=[], loglikelihoods=[], num_params=[], aics=[])
            tests = ContinuousTests(
                names=["Test 1", "Test 2", "Test 3", "Test 4"],
                ll_ratios=[constants.BMDS_BLANK_VALUE] * 4,
                dfs=[constants.BMDS_BLANK_VALUE] * 4,
                p_values=[constants.BMDS_BLANK_VALUE] * 4,
            )
        else:
            fit = model.results.fit.model_copy(update={"bmd_dist": draws})
            gof = (
                ContinuousGof.from_cpp(cpp_result.gof, summary.bmd, model.dataset.doses)
                if cpp_result is not None
                else model.results.gof
            )
            deviance = model.results.deviance
            tests = model.results.tests

        p_value = (
            self.model_p_values[index]
            if index < len(self.model_p_values)
            else constants.BMDS_BLANK_VALUE
        )
        waic = (
            self.model_waics[index] if index < len(self.model_waics) else constants.BMDS_BLANK_VALUE
        )
        loglikelihood = (
            self.model_loglikelihoods[index]
            if index < len(self.model_loglikelihoods)
            else constants.BMDS_BLANK_VALUE
        )
        fit = fit.model_copy(update={"loglikelihood": loglikelihood})

        extra_values = [summary.bmd] if summary.bmd >= 0 else []
        dr_x = model.dataset.dose_linspace(extra_values=extra_values)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*invalid value encountered.*")
            warnings.filterwarnings("ignore", message=".*divide by zero encountered.*")
            dr_y = clean_array(model.dr_curve(dr_x, parameters.values))

        xs = np.asarray([summary.bmdl, summary.bmd, summary.bmdu], dtype=float)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", message=".*invalid value encountered.*")
            warnings.filterwarnings("ignore", message=".*divide by zero encountered.*")
            critical_ys = clean_array(model.dr_curve(xs, parameters.values))
        critical_ys[critical_ys <= 0] = constants.BMDS_BLANK_VALUE

        plotting = ContinuousPlotting(
            dr_x=dr_x,
            dr_y=dr_y,
            bmdl_y=float(critical_ys[0]),
            bmd_y=float(critical_ys[1]),
            bmdu_y=float(critical_ys[2]),
        )

        model.results = ContinuousResult(
            bmdl=summary.bmdl,
            bmd=summary.bmd,
            bmdu=summary.bmdu,
            has_completed=True,
            fit=fit,
            gof=gof,
            parameters=parameters,
            deviance=deviance,
            tests=tests,
            plotting=plotting,
            summary_p_value=p_value,
            summary_waic=waic,
        )
