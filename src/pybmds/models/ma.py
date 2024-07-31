from pydantic import Field

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
        self.structs = DichotomousModelAverage(
            self.session.dataset, self.models, self.session.ma_weights
        )
        self.structs.execute()
        return DichotomousModelAverageResult.from_cpp(
            self.structs, [model.results for model in self.models]
        )

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
