import numpy as np
from pydantic import Field

from .. import bmdscore, plotting
from ..constants import ZEROISH, NestedDichotomousModel, NestedDichotomousModelChoices, PriorClass
from ..datasets import NestedDichotomousDataset
from ..types.nested_dichotomous import (
    NestedDichotomousAnalysis,
    NestedDichotomousModelSettings,
    NestedDichotomousResult,
)
from ..types.priors import ModelPriors, get_nested_dichotomous_prior
from ..utils import multi_lstrip
from .base import BmdModel, BmdModelSchema, InputModelSettings


class BmdModelNestedDichotomous(BmdModel):
    bmd_model_class: NestedDichotomousModel
    model_class: bmdscore.nested_model

    def name(self) -> str:
        return (
            self.settings.name
            or f"{super().name()} ({self.settings.litter_specific_covariate.text}{self.settings.intralitter_correlation.text})"
        )

    def get_model_settings(
        self, dataset: NestedDichotomousDataset, settings: InputModelSettings
    ) -> NestedDichotomousModelSettings:
        if settings is None:
            model_settings = NestedDichotomousModelSettings()
        elif isinstance(settings, NestedDichotomousModelSettings):
            model_settings = settings
        else:
            model_settings = NestedDichotomousModelSettings.model_validate(settings)

        # get default values, may require further model customization
        if not isinstance(model_settings.priors, ModelPriors):
            prior_class = (
                model_settings.priors
                if isinstance(model_settings.priors, PriorClass)
                else self.get_default_prior_class()
            )
            model_settings.priors = get_nested_dichotomous_prior(
                self.bmd_model_class, prior_class=prior_class
            )

        return model_settings

    def get_default_prior_class(self) -> PriorClass:
        return PriorClass.frequentist_restricted

    def to_cpp(self) -> NestedDichotomousAnalysis:
        return NestedDichotomousAnalysis.create(
            model_class=self.model_class,
            nparms=len(self.get_param_names()),
            dataset=self.dataset,
            settings=self.settings,
        )

    def execute(self) -> NestedDichotomousResult:
        self.structs = self.to_cpp()
        self.structs.execute()
        self.results = NestedDichotomousResult.from_model(self)
        return self.results

    def _plot_bmr_lines(self, ax, axlines: bool):
        plotting.add_bmr_lines(
            ax,
            self.results.bmd,
            self.results.plotting.bmd_y,
            self.results.bmdl,
            self.results.bmdu,
            axlines=axlines,
        )

    def serialize(self) -> "BmdModelNestedDichotomousSchema":
        return BmdModelNestedDichotomousSchema(
            name=self.name(),
            model_class=self.bmd_model_class,
            settings=self.settings,
            results=self.results,
        )

    def get_gof_pvalue(self): ...

    def get_priors_list(self) -> list[list]:
        return self.settings.priors.priors_list(nphi=self.dataset.num_dose_groups)

    def model_settings_text(self) -> str:
        input_tbl = self.settings.tbl()
        return multi_lstrip(
            f"""
        Input Summary:
        {input_tbl}
        """
        )


class BmdModelNestedDichotomousSchema(BmdModelSchema):
    name: str
    bmd_model_class: NestedDichotomousModel = Field(alias="model_class")
    settings: NestedDichotomousModelSettings
    results: NestedDichotomousResult | None

    def deserialize(self, dataset: NestedDichotomousDataset) -> BmdModelNestedDichotomous:
        Model = bmd_model_map[self.bmd_model_class.id]
        model = Model(dataset=dataset, settings=self.settings)
        model.results = self.results
        return model


class NestedLogistic(BmdModelNestedDichotomous):
    bmd_model_class = NestedDichotomousModelChoices.logistic.value
    model_class = bmdscore.nested_model.nlogistic

    def get_param_names(self) -> list[str]:
        return ["a", "b", "theta1", "theta2", "rho"] + [
            f"phi{i}" for i in range(1, self.dataset.num_dose_groups + 1)
        ]

    def dr_curve(self, doses: np.ndarray, params: dict, fixed_lsc: float) -> np.ndarray:
        a = params["a"]
        b = params["b"]
        theta1 = params["theta1"]
        theta2 = params["theta2"]
        rho = params["rho"]
        d = doses.copy()
        d[d < ZEROISH] = ZEROISH
        return (
            a
            + theta1 * fixed_lsc
            + (1 - a - theta1 * fixed_lsc)
            / (1 + np.exp(-1 * b - theta2 * fixed_lsc - rho * np.log(d)))
        )

    def get_default_prior_class(self) -> PriorClass:
        return PriorClass.frequentist_restricted


class Nctr(BmdModelNestedDichotomous):
    bmd_model_class = NestedDichotomousModelChoices.nctr.value
    model_class = bmdscore.nested_model.nctr

    def get_param_names(self) -> list[str]:
        return ["a", "b", "theta1", "theta2", "rho"] + [
            f"phi{i}" for i in range(1, self.dataset.num_dose_groups + 1)
        ]

    def dr_curve(self, doses: np.ndarray, params: dict, fixed_lsc: float) -> np.ndarray:
        raise NotImplementedError("TODO - update formula")

    def get_default_prior_class(self) -> PriorClass:
        return PriorClass.frequentist_restricted

    def get_model_settings(
        self, dataset: NestedDichotomousDataset, settings: InputModelSettings
    ) -> NestedDichotomousModelSettings:
        model_settings = super().get_model_settings(dataset, settings)

        smax = max(1, max(dataset.litter_ns))
        smin = max(1, min(dataset.litter_ns))
        model_settings.priors.update("theta1", min_value=-1.0 / smax, max_value=-1.0 / smin)
        model_settings.priors.update("theta2", min_value=-1.0 / smax, max_value=-1.0 / smin)

        return model_settings


bmd_model_map = {
    NestedDichotomousModelChoices.logistic.value.id: NestedLogistic,
    NestedDichotomousModelChoices.nctr.value.id: Nctr,
}
