import numpy as np
from pydantic import Field

from ..constants import ContinuousModel, ContinuousModelChoices, DistType, PriorClass
from ..datasets.continuous import ContinuousDatasets
from ..exceptions import ConfigurationException
from ..types.continuous import ContinuousAnalysis, ContinuousModelSettings, ContinuousResult
from ..types.priors import ModelPriors, get_continuous_prior
from .base import BmdModel, BmdModelSchema, InputModelSettings


class BmdModelContinuous(BmdModel):
    bmd_model_class: ContinuousModel

    def get_model_settings(
        self,
        dataset: ContinuousDatasets,
        settings: InputModelSettings,
    ) -> ContinuousModelSettings:
        if settings is None:
            model_settings = ContinuousModelSettings()
        elif isinstance(settings, ContinuousModelSettings):
            model_settings = settings
        else:
            model_settings = ContinuousModelSettings.model_validate(settings)

        # only estimate direction if unspecified
        if model_settings.is_increasing is None:
            model_settings.is_increasing = dataset.is_increasing

        # set preferred variance model if unspecified
        if settings is None or (isinstance(settings, dict) and "disttype" not in settings):
            model_settings.disttype = self.set_constant_variance_value()

        # get default values, may require further model customization
        if not isinstance(model_settings.priors, ModelPriors):
            prior_class = (
                model_settings.priors
                if isinstance(model_settings.priors, PriorClass)
                else self.get_default_prior_class()
            )
            model_settings.priors = get_continuous_prior(
                self.bmd_model_class, prior_class=prior_class
            )

        return model_settings

    def set_constant_variance_value(self) -> DistType:
        # set modeled variance if p-test 2 < 0.05 or not calculated, otherwise constant
        anova = self.dataset.anova()
        return (
            DistType.normal_ncv if (anova is None or anova.test2.TEST < 0.05) else DistType.normal
        )

    def _build_inputs(self) -> ContinuousAnalysis:
        return ContinuousAnalysis(
            model=self.bmd_model_class,
            dataset=self.dataset,
            priors=self.settings.priors,
            BMD_type=self.settings.bmr_type,
            BMR=self.settings.bmr,
            alpha=self.settings.alpha,
            is_increasing=self.settings.is_increasing,
            tail_prob=self.settings.tail_prob,
            disttype=self.settings.disttype,
            samples=self.settings.samples,
            burnin=self.settings.burnin,
            degree=self.settings.degree,
            count_all_parameters_on_boundary=self.settings.count_all_parameters_on_boundary,
        )

    def execute(self) -> ContinuousResult:
        inputs = self._build_inputs()
        structs = inputs.to_cpp()
        self.structs = structs
        self.structs.execute()
        self.results = ContinuousResult.from_model(self)
        return self.results

    def get_default_model_degree(self, dataset) -> int:
        return 0

    def get_default_prior_class(self) -> PriorClass:
        return PriorClass.frequentist_restricted

    def serialize(self) -> "BmdModelContinuousSchema":
        return BmdModelContinuousSchema(
            name=self.name(),
            model_class=self.bmd_model_class,
            settings=self.settings,
            results=self.results,
        )

    def get_param_names(self) -> list[str]:
        names = list(self.bmd_model_class.params)
        names.extend(self.get_variance_param_names())
        return names

    def get_variance_param_names(self):
        if self.settings.disttype == DistType.normal_ncv:
            return list(self.bmd_model_class.variance_params)
        else:
            return [self.bmd_model_class.variance_params[1]]

    def get_gof_pvalue(self) -> float:
        return self.results.tests.p_values[3]

    def get_priors_list(self) -> list[list]:
        degree = self.settings.degree if self.degree_required else None
        return self.settings.priors.priors_list(degree=degree, dist_type=self.settings.disttype)


class BmdModelContinuousSchema(BmdModelSchema):
    name: str
    bmds_model_class: ContinuousModel = Field(alias="model_class")
    settings: ContinuousModelSettings
    results: ContinuousResult | None = None

    def deserialize(self, dataset: ContinuousDatasets) -> BmdModelContinuous:
        Model = get_model_class(self)
        model = Model(dataset=dataset, settings=self.settings)
        model.results = self.results
        return model


class Power(BmdModelContinuous):
    bmd_model_class = ContinuousModelChoices.power.value

    def get_model_settings(
        self, dataset: ContinuousDatasets, settings: InputModelSettings
    ) -> ContinuousModelSettings:
        model_settings = super().get_model_settings(dataset, settings)

        if model_settings.disttype == DistType.log_normal:
            raise ConfigurationException("Power model cannot run with lognormal distribution")

        if model_settings.priors.prior_class in [
            PriorClass.frequentist_unrestricted,
            PriorClass.frequentist_restricted,
        ]:
            is_cv = model_settings.disttype in [DistType.normal, DistType.log_normal]
            v = model_settings.priors.get_prior("v")
            v.min_value = -100 if is_cv else -10_000
            v.max_value = 100 if is_cv else 10_000

        return model_settings

    def dr_curve(self, doses, params) -> np.ndarray:
        g = params[0]
        v = params[1]
        n = params[2]
        return g + v * doses**n


class Hill(BmdModelContinuous):
    bmd_model_class = ContinuousModelChoices.hill.value

    def dr_curve(self, doses, params) -> np.ndarray:
        g = params[0]
        v = params[1]
        k = params[2]
        n = params[3]
        return g + v * doses**n / (k**n + doses**n)

    def get_model_settings(
        self, dataset: ContinuousDatasets, settings: InputModelSettings
    ) -> ContinuousModelSettings:
        model_settings = super().get_model_settings(dataset, settings)

        if model_settings.disttype == DistType.log_normal:
            raise ConfigurationException("Hill model cannot run with lognormal distribution")

        return model_settings


class Polynomial(BmdModelContinuous):
    bmd_model_class = ContinuousModelChoices.polynomial.value
    degree_required: bool = True

    def name(self) -> str:
        return self.settings.name or f"Polynomial {self.settings.degree}"

    def get_model_settings(
        self, dataset: ContinuousDatasets, settings: InputModelSettings
    ) -> ContinuousModelSettings:
        model_settings = super().get_model_settings(dataset, settings)

        if model_settings.disttype == DistType.log_normal:
            raise ConfigurationException("Polynomial model cannot run with lognormal distribution")

        if model_settings.degree < 1:
            model_settings.degree = self.get_default_model_degree(dataset)

        is_freq = model_settings.priors.prior_class in [
            PriorClass.frequentist_restricted,
            PriorClass.frequentist_unrestricted,
        ]
        if is_freq:
            g = model_settings.priors.get_prior("g")
            b1 = model_settings.priors.get_prior("b1")
            bN = model_settings.priors.get_prior("bN")
            is_cv = model_settings.disttype in [DistType.normal, DistType.log_normal]
            # update mins
            g.min_value = -1e6 if is_cv else 0
            b1.min_value = -1e6 if is_cv else -18
            # update maxes
            g.max_value = 1e6 if is_cv else 1_000
            b1.max_value = 1e6 if is_cv else 18
            # for restricted, betas in one direction
            if model_settings.priors.prior_class is PriorClass.frequentist_restricted:
                attr = "min_value" if model_settings.is_increasing else "max_value"
                setattr(b1, attr, 0)
                setattr(bN, attr, 0)

        return model_settings

    def get_default_model_degree(self, dataset) -> int:
        return 2

    def dr_curve(self, doses, params) -> np.ndarray:
        val = doses * 0.0 + params[0]
        for i in range(1, self.settings.degree + 1):
            val += params[i] * doses**i
        return val

    def get_param_names(self) -> list[str]:
        names = [f"b{i}" for i in range(self.settings.degree + 1)]
        names.extend(self.get_variance_param_names())
        names[0] = "g"
        return names


class Linear(Polynomial):
    def name(self) -> str:
        return self.settings.name or "Linear"

    def get_default_model_degree(self, dataset) -> int:
        return 1

    def get_default_prior_class(self) -> PriorClass:
        return PriorClass.frequentist_unrestricted

    def get_model_settings(
        self, dataset: ContinuousDatasets, settings: InputModelSettings
    ) -> ContinuousModelSettings:
        model_settings = super().get_model_settings(dataset, settings)

        if model_settings.disttype == DistType.log_normal:
            raise ConfigurationException("Linear model cannot run with lognormal distribution")

        if model_settings.degree != 1:
            raise ConfigurationException("Linear model must have degree of 1")

        return model_settings


class ExponentialM3(BmdModelContinuous):
    bmd_model_class = ContinuousModelChoices.exp_m3.value

    def get_model_settings(
        self, dataset: ContinuousDatasets, settings: InputModelSettings
    ) -> ContinuousModelSettings:
        model_settings = super().get_model_settings(dataset, settings)

        if model_settings.priors.prior_class is PriorClass.frequentist_restricted:
            attr = "min_value" if model_settings.is_increasing else "max_value"
            setattr(model_settings.priors.get_prior("c"), attr, 0.0)

        return model_settings

    def dr_curve(self, doses, params) -> np.ndarray:
        a = params[0]
        b = params[1]
        d = params[2]
        sign = 1.0 if self.settings.is_increasing else -1.0
        return a * np.exp(sign * ((b * doses) ** d))


class ExponentialM5(BmdModelContinuous):
    bmd_model_class = ContinuousModelChoices.exp_m5.value

    def get_model_settings(
        self, dataset: ContinuousDatasets, settings: InputModelSettings
    ) -> ContinuousModelSettings:
        model_settings = super().get_model_settings(dataset, settings)

        if model_settings.priors.prior_class is PriorClass.frequentist_restricted:
            attr = "min_value" if model_settings.is_increasing else "max_value"
            setattr(model_settings.priors.get_prior("c"), attr, 0.0)

        return model_settings

    def dr_curve(self, doses, params) -> np.ndarray:
        a = params[0]
        b = params[1]
        c = params[2]
        d = params[3]
        return a * (c - (c - 1.0) * np.exp(-1.0 * np.power(b * doses, d)))


_bmd_model_map = {
    ContinuousModelChoices.power.value.id: Power,
    ContinuousModelChoices.hill.value.id: Hill,
    ContinuousModelChoices.exp_m3.value.id: ExponentialM3,
    ContinuousModelChoices.exp_m5.value.id: ExponentialM5,
}


def get_model_class(data: BmdModelContinuousSchema) -> type[BmdModelContinuous]:
    """Get continuous model class given the schema

    Generally this is a dictionary lookup; however there is a special case for
    Linear/Polynomial because they have the same model class enum in C++, but different
    python classes. Thus we specify by the degree as well.
    """
    if data.bmds_model_class.id == ContinuousModelChoices.polynomial.value.id:
        return Linear if data.settings.degree == 1 else Polynomial
    else:
        return _bmd_model_map[data.bmds_model_class.id]
