import numpy as np
from pydantic import Field

from ..constants import BMDS_BLANK_VALUE, PriorClass
from ..types.cma import ContinuousModelAverage, ContinuousModelAverageResult
from ..types.continuous import ContinuousModelSettings
from .base import BmdModelAveraging, BmdModelAveragingSchema, InputModelSettings


class BmdModelAveragingContinuous(BmdModelAveraging):
    def get_model_settings(self, settings: InputModelSettings) -> ContinuousModelSettings:
        if settings is None:
            return ContinuousModelSettings()
        elif isinstance(settings, ContinuousModelSettings):
            return settings
        else:
            return ContinuousModelSettings.model_validate(settings)

    def execute(self) -> ContinuousModelAverageResult:
        model_indexes = [self.session.models.index(model) for model in self.models]
        model_weights = np.asarray(self.session.ma_weights, dtype=float)
        if model_weights.size != len(self.models):
            model_weights = model_weights[model_indexes]
            model_weights = model_weights / model_weights.sum()

        self.structs = ContinuousModelAverage(
            self.session.dataset, self.models, model_weights, self.session.weight_option
        )
        avg = self.structs.average  # python_continuousMA_analysis
        n = int(avg.nmodels)

        # NOTE: your current binding exposes loud_dist_type as "disttype"
        if len(avg.models) != n:
            raise ValueError(f"avg.models length {len(avg.models)} != nmodels {n}")
        if not hasattr(avg, "disttype"):
            raise AttributeError("avg has no 'disttype' (bound loud_dist_type)")
        if len(avg.disttype) != n:
            raise ValueError(f"avg.disttype length {len(avg.disttype)} != nmodels {n}")
        self.structs.execute()
        results = ContinuousModelAverageResult.from_cpp(self.structs)
        for idx, model in enumerate(self.models):
            if model.settings.priors.prior_class is PriorClass.bayesian_loud:
                results.sync_model_result(model, idx, self.structs.result.models[idx])
        extra_values = [results.bmd] if results.bmd >= 0 else []
        dr_x = self.session.dataset.dose_linspace(extra_values=extra_values)
        values = np.asarray(
            [
                np.interp(dr_x, model.results.plotting.dr_x, model.results.plotting.dr_y)
                for model in self.models
            ],
            dtype=float,
        )
        posteriors = np.asarray(results.posteriors, dtype=float)
        valid_posteriors = (
            np.isfinite(posteriors) & (posteriors >= 0) & (posteriors != BMDS_BLANK_VALUE)
        )
        if values.size and valid_posteriors.any() and posteriors[valid_posteriors].sum() > 0:
            weights = posteriors.copy()
            weights[~valid_posteriors] = 0.0
            weights = weights / weights.sum()
            dr_y = weights @ values
            bmds = np.asarray([results.bmdl, results.bmd, results.bmdu], dtype=float)
            bmds_ys = np.interp(bmds, dr_x, dr_y)
            results = results.model_copy(
                update={
                    "dr_x": dr_x,
                    "dr_y": dr_y,
                    "bmdl_y": float(bmds_ys[0]),
                    "bmd_y": float(bmds_ys[1]),
                    "bmdu_y": float(bmds_ys[2]),
                }
            )
        return results

    def serialize(self, session) -> "BmdModelAveragingContinuousSchema":
        model_indexes = [session.models.index(model) for model in self.models]
        return BmdModelAveragingContinuousSchema(
            settings=self.settings, model_indexes=model_indexes, results=self.results
        )


class BmdModelAveragingContinuousSchema(BmdModelAveragingSchema):
    settings: ContinuousModelSettings
    results: ContinuousModelAverageResult
    bmds_model_indexes: list[int] = Field(alias="model_indexes")

    def deserialize(self, session) -> BmdModelAveragingContinuous:
        models = [session.models[idx] for idx in self.bmds_model_indexes]
        ma = BmdModelAveragingContinuous(session=session, models=models, settings=self.settings)
        ma.results = self.results
        return ma
