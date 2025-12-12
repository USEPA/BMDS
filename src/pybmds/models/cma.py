from pydantic import Field

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
        self.structs = ContinuousModelAverage(
            self.session.dataset, self.models, self.session.ma_weights
        )
        self.structs.execute()
        return ContinuousModelAverageResult.from_cpp(
            self.structs, [model.results for model in self.models]
        )

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
