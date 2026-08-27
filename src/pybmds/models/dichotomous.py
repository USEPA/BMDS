from typing import Self

import numpy as np
from pydantic import Field
from scipy.stats import gamma, norm

from ..constants import DichotomousModel, DichotomousModelChoices, PriorClass
from ..datasets import DichotomousDataset
from ..types.dichotomous import DichotomousAnalysis, DichotomousModelSettings, DichotomousResult
from ..types.priors import ModelPriors, get_dichotomous_prior, multistage_cancer_prior
from ..utils import multi_lstrip
from .base import BmdModel, BmdModelSchema, InputModelSettings


class BmdModelDichotomous(BmdModel):
    bmd_model_class: DichotomousModel

    def results_text(self) -> str:
        if self.results is None:  # pragma: no cover
            raise ValueError("Cannot render text if results are unavailable")

        session = getattr(self, "session", None)
        if self.settings.priors.prior_class is PriorClass.bayesian_loud and session is not None:
            summary = session.get_model_average_summary_for_model(self)
            if summary is not None:
                return multi_lstrip(
                    f"""
                Modeling Summary:
                {self.results.tbl()}

                Model Parameters:
                {self.results.parameters.tbl()}

                Goodness of Fit:
                {self.results.gof.tbl(self.dataset)}

                LOUD Model-Average Weights:
                Prior Weight: {summary.prior}
                Posterior Weight: {summary.posterior}

                Model-specific BMD values shown above are taken from the LOUD
                model averaging result for this model.
                """
                )

        return super().results_text()

    def get_model_settings(
        self, dataset: DichotomousDataset, settings: InputModelSettings
    ) -> DichotomousModelSettings:
        if settings is None:
            model_settings = DichotomousModelSettings()
        elif isinstance(settings, DichotomousModelSettings):
            model_settings = settings
        else:
            model_settings = DichotomousModelSettings.model_validate(settings)

        # get default values, may require further model customization
        if not isinstance(model_settings.priors, ModelPriors):
            prior_class = (
                model_settings.priors
                if isinstance(model_settings.priors, PriorClass)
                else self.get_default_prior_class()
            )
            model_settings.priors = get_dichotomous_prior(
                self.bmd_model_class, prior_class=prior_class
            )
        if model_settings.priors.prior_class is PriorClass.bayesian_loud:
            model_settings = DichotomousModelSettings.model_validate(
                {
                    field: getattr(model_settings, field)
                    for field in DichotomousModelSettings.model_fields
                }
            )

        return model_settings

    def _build_inputs(self) -> DichotomousAnalysis:
        return DichotomousAnalysis(
            model=self.bmd_model_class,
            dataset=self.dataset,
            priors=self.settings.priors,
            BMD_type=self.settings.bmr_type,
            BMR=self.settings.bmr,
            alpha=self.settings.alpha,
            degree=self.settings.degree,
            samples=self.settings.samples,
            burnin=self.settings.burnin,
            n_chains=self.settings.n_chains,
            seed=self.settings.seed,
            count_all_parameters_on_boundary=self.settings.count_all_parameters_on_boundary,
        )

    def execute(self, slope_factor: bool = False) -> DichotomousResult:
        """Execute analysis using bmdscore

        Args:
            slope_factor (bool, optional; default False): If True, calculates slope factor

        Returns:
            DichotomousResult: _description_
        """
        if self.settings.priors.prior_class is PriorClass.bayesian_loud and (
            self.session is None or self.session.model_average is None
        ):
            return self._execute_standalone_loud()

        inputs = self._build_inputs()
        structs = inputs.to_cpp()
        self.structs = structs
        if (
            self.settings.priors.prior_class is PriorClass.bayesian_loud
            and getattr(self, "session", None) is not None
            and self.session.model_average is not None
        ):
            self.results = None
            return self.results
        self.structs.execute()
        if slope_factor:
            bmr = self.structs.analysis.BMR
            self.structs.result.bmdsRes.setSlopeFactor(bmr)
        self.results = DichotomousResult.from_model(self)
        return self.results

    def get_default_model_degree(self, dataset) -> int:
        return 2

    def get_default_prior_class(self) -> PriorClass:
        return PriorClass.frequentist_restricted

    def get_param_names(self) -> list[str]:
        names = list(self.bmd_model_class.params)
        return names

    def serialize(self) -> "BmdModelDichotomousSchema":
        return BmdModelDichotomousSchema(
            name=self.name(),
            model_class=self.bmd_model_class,
            settings=self.settings,
            results=self.results,
        )

    def get_gof_pvalue(self) -> float:
        return self.results.gof.p_value

    def get_priors_list(self) -> list[list]:
        degree = self.settings.degree if self.degree_required else None
        return self.settings.priors.priors_list(degree=degree)


class BmdModelDichotomousSchema(BmdModelSchema):
    name: str
    bmds_model_class: DichotomousModel = Field(alias="model_class")
    settings: DichotomousModelSettings
    results: DichotomousResult | None = None

    def deserialize(self, dataset: DichotomousDataset) -> BmdModelDichotomous:
        Model = bmd_model_map[self.bmds_model_class.id]
        model = Model(dataset=dataset, settings=self.settings)
        model.results = self.results
        return model


class Logistic(BmdModelDichotomous):
    bmd_model_class = DichotomousModelChoices.logistic.value

    def dr_curve(self, doses, params) -> np.ndarray:
        a = params[0]
        b = params[1]
        return 1 / (1 + np.exp(-a - b * doses))

    def get_default_prior_class(self) -> PriorClass:
        return PriorClass.frequentist_unrestricted


class LogLogistic(BmdModelDichotomous):
    bmd_model_class = DichotomousModelChoices.loglogistic.value

    def dr_curve(self, doses, params) -> np.ndarray:
        g = params[0]
        a = params[1]
        b = params[2]
        return g + (1 - g) * (1 / (1 + np.exp(-a - b * np.log(doses))))


class Probit(BmdModelDichotomous):
    bmd_model_class = DichotomousModelChoices.probit.value

    def dr_curve(self, doses, params) -> np.ndarray:
        a = params[0]
        b = params[1]
        return norm.cdf(a + b * doses)

    def get_default_prior_class(self) -> PriorClass:
        return PriorClass.frequentist_unrestricted


class LogProbit(BmdModelDichotomous):
    bmd_model_class = DichotomousModelChoices.logprobit.value

    def dr_curve(self, doses, params) -> np.ndarray:
        g = params[0]
        a = params[1]
        b = params[2]
        return g + (1 - g) * norm.cdf(a + b * np.log(doses))

    def get_default_prior_class(self) -> PriorClass:
        return PriorClass.frequentist_unrestricted


class Gamma(BmdModelDichotomous):
    bmd_model_class = DichotomousModelChoices.gamma.value

    def dr_curve(self, doses, params) -> np.ndarray:
        g = params[0]
        a = params[1]
        b = params[2]
        return g + (1 - g) * gamma.cdf(b * doses, a)


class QuantalLinear(BmdModelDichotomous):
    bmd_model_class = DichotomousModelChoices.qlinear.value

    def dr_curve(self, doses, params) -> np.ndarray:
        g = params[0]
        b = params[1]
        return g + (1 - g) * (1 - np.exp(-b * doses))

    def get_default_prior_class(self) -> PriorClass:
        return PriorClass.frequentist_unrestricted


class Weibull(BmdModelDichotomous):
    bmd_model_class = DichotomousModelChoices.weibull.value

    def dr_curve(self, doses, params) -> np.ndarray:
        g = params[0]
        a = params[1]
        b = params[2]
        return g + (1 - g) * (1 - np.exp(-b * doses**a))


class DichotomousHill(BmdModelDichotomous):
    bmd_model_class = DichotomousModelChoices.hill.value

    def dr_curve(self, doses, params) -> np.ndarray:
        g = params[0]
        v = params[1]
        a = params[2]
        b = params[3]
        return g + (1 - g) * v * (1 / (1 + np.exp(-a - b * np.log(doses))))


class Multistage(BmdModelDichotomous):
    bmd_model_class = DichotomousModelChoices.multistage.value
    degree_required: bool = True

    def get_model_settings(
        self, dataset: DichotomousDataset, settings: InputModelSettings
    ) -> DichotomousModelSettings:
        model_settings = super().get_model_settings(dataset, settings)

        if model_settings.degree < 1:
            model_settings.degree = self.get_default_model_degree(dataset)

        return model_settings

    def name(self) -> str:
        return self.settings.name or f"Multistage {self.settings.degree}"

    def dr_curve(self, doses, params) -> np.ndarray:
        g = params[0]
        val = doses * 0
        for i in range(1, len(params)):
            val += params[i] * doses**i
        return g + (1 - g) * (1 - np.exp(-1.0 * val))

    def get_param_names(self) -> list[str]:
        names = [f"b{i}" for i in range(self.settings.degree + 1)]
        names[0] = "g"
        return names


class MultistageCancer(Multistage):
    def get_model_settings(
        self, dataset: DichotomousDataset, settings: InputModelSettings
    ) -> DichotomousModelSettings:
        override_default_prior = settings is None or (
            isinstance(settings, dict) and "priors" not in settings
        )
        model_settings = super().get_model_settings(dataset, settings)
        if override_default_prior:
            model_settings.priors = self.custom_prior()
        return model_settings

    def custom_prior(self) -> ModelPriors:
        return multistage_cancer_prior()

    @classmethod
    def deserialize(cls, dataset: DichotomousDataset, obj: BmdModelDichotomousSchema) -> Self:
        model = cls(dataset=dataset, settings=obj.settings)
        model.results = obj.results
        return model


bmd_model_map = {
    DichotomousModelChoices.hill.value.id: DichotomousHill,
    DichotomousModelChoices.gamma.value.id: Gamma,
    DichotomousModelChoices.logistic.value.id: Logistic,
    DichotomousModelChoices.loglogistic.value.id: LogLogistic,
    DichotomousModelChoices.logprobit.value.id: LogProbit,
    DichotomousModelChoices.multistage.value.id: Multistage,
    DichotomousModelChoices.probit.value.id: Probit,
    DichotomousModelChoices.qlinear.value.id: QuantalLinear,
    DichotomousModelChoices.weibull.value.id: Weibull,
}
