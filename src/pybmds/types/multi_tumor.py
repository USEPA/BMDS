from typing import NamedTuple, Self

from pydantic import BaseModel, ConfigDict, Field

from .. import bmdscore
from ..models.dichotomous import BmdModelDichotomousSchema
from ..utils import multi_lstrip, pretty_table
from .common import inspect_cpp_obj
from .dichotomous import (
    DichotomousAnalysisCPPStructs,
    DichotomousResult,
    DichotomousRiskType,
)


class MultitumorAnalysis(NamedTuple):
    analysis: bmdscore.python_multitumor_analysis
    result: bmdscore.python_multitumor_result

    def execute(self):
        bmdscore.pythonBMDSMultitumor(self.analysis, self.result)

    def __str__(self) -> str:
        lines = []
        inspect_cpp_obj(lines, self.analysis, depth=0)
        inspect_cpp_obj(lines, self.result, depth=0)
        return "\n".join(lines)


class MultitumorSettings(BaseModel):
    degrees: list[int]
    bmr: float = Field(default=0.1, gt=0, lt=1)
    alpha: float = Field(default=0.05, gt=0, lt=1)
    bmr_type: DichotomousRiskType = DichotomousRiskType.ExtraRisk

    model_config = ConfigDict(extra="forbid")


class MultitumorResult(BaseModel):
    bmd: float
    bmdl: float
    bmdu: float
    ll: float
    ll_constant: float
    models: list[list[BmdModelDichotomousSchema]]  # all degrees for all datasets
    selected_model_indexes: list[int | None]
    slope_factor: float
    valid_result: bool

    @classmethod
    def from_model(cls, model) -> Self:
        result: bmdscore.python_multitumor_result = model.structs.result
        i_models = []
        for i, models in enumerate(model.models):
            j_models = []
            i_models.append(j_models)
            for j, m in enumerate(models):
                m.structs = DichotomousAnalysisCPPStructs(
                    analysis=model.structs.analysis.models[i][j],
                    result=model.structs.result.models[i][j],
                )
                m.results = DichotomousResult.from_model(m)
                j_models.append(m.serialize())

        return cls(
            bmd=result.BMD,
            bmdl=result.BMDL,
            bmdu=result.BMDU,
            ll=result.combined_LL,
            ll_constant=result.combined_LL_const,
            models=i_models,
            selected_model_indexes=[idx if idx >= 0 else None for idx in result.selectedModelIndex],
            slope_factor=result.slopeFactor,
            valid_result=result.validResult,
        )

    def text(self, datasets, models) -> str:
        texts = []
        for i, dataset in enumerate(datasets):
            model_idx = self.selected_model_indexes[i]
            texts.append("\n" + dataset._get_dataset_name() + "\n" + "═" * 80)
            texts.append("\n" + dataset.tbl() + "\n")
            if model_idx is None:
                texts.append("No model selected.")
            else:
                texts.append(models[i][model_idx].text())
        fitted = "\n".join(texts)

        return multi_lstrip(
            f"""
        Modeling Summary:
        {self.tbl()}

        {fitted}
        """
        )

    def tbl(self) -> str:
        data = [
            ["BMD", self.bmd],
            ["BMDL", self.bmdl],
            ["BMDU", self.bmdu],
            ["Slope Factor", self.slope_factor],
            ["Combined Log-likelihood", self.ll],
            ["Combined Log-likelihood Constant", self.ll_constant],
        ]
        return pretty_table(data, "")
