import numpy as np
from pydantic import Field

from ..constants import PriorClass
from ..types.dichotomous import DichotomousModelSettings
from ..types.ma import DichotomousModelAverage, DichotomousModelAverageResult
from .base import BmdModelAveraging, BmdModelAveragingSchema, InputModelSettings


class BmdModelAveragingDichotomous(BmdModelAveraging):
    def get_model_settings(self, settings: InputModelSettings) -> DichotomousModelSettings:
        if settings is None:
            return DichotomousModelSettings()
        elif isinstance(settings, DichotomousModelSettings):
            return settings
        else:
            return DichotomousModelSettings.model_validate(settings)

    def execute(self) -> DichotomousModelAverageResult:
        is_loud = self.settings.priors.prior_class is PriorClass.bayesian_loud
        model_weights = self.session.ma_weights
        if is_loud:
            model_indexes = [self.session.models.index(model) for model in self.models]
            model_weights = np.asarray(model_weights, dtype=float)
            if model_weights.size == len(self.session.models):
                model_weights = model_weights[model_indexes]
                model_weights = model_weights / model_weights.sum()
            elif model_weights.size != len(self.models):
                model_weights = model_weights[model_indexes]
                model_weights = model_weights / model_weights.sum()

        self.structs = DichotomousModelAverage(
            self.session.dataset,
            self.models,
            model_weights,
            self.settings.priors.prior_class,
            self.session.weight_option,
        )
        self.structs.execute()
        results = DichotomousModelAverageResult.from_cpp(
            self.structs, self.models if is_loud else [model.results for model in self.models]
        )
        if is_loud:
            for idx, model in enumerate(self.models):
                results.sync_model_result(model, idx, self.structs.result.models[idx])
        return results

    def serialize(self, session) -> "BmdModelAveragingDichotomousSchema":
        model_indexes = [session.models.index(model) for model in self.models]
        return BmdModelAveragingDichotomousSchema(
            settings=self.settings, model_indexes=model_indexes, results=self.results
        )


class BmdModelAveragingDichotomousSchema(BmdModelAveragingSchema):
    settings: DichotomousModelSettings
    results: DichotomousModelAverageResult
    bmds_model_indexes: list[int] = Field(alias="model_indexes")

    def deserialize(self, session) -> BmdModelAveragingDichotomous:
        models = [session.models[idx] for idx in self.bmds_model_indexes]
        ma = BmdModelAveragingDichotomous(session=session, models=models, settings=self.settings)
        ma.results = self.results
        return ma
