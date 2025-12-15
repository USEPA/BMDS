from enum import Enum, IntEnum, StrEnum

from pydantic import BaseModel, Field

from .bmdscore import cont_model, dich_model, nested_model, loud_cont_model

ZEROISH = 1e-8
BMDS_BLANK_VALUE = -9999
N_BMD_DIST = 100
NUM_PRIOR_COLS = 5
MAXIMUM_POLYNOMIAL_ORDER = 8
BOOL_YES_NO = {True: "yes", False: "no"}


class ModelClass(StrEnum):
    # Types of modeling sessions
    DICHOTOMOUS = "D"
    CONTINUOUS = "C"
    NESTED_DICHOTOMOUS = "ND"
    MULTI_TUMOR = "MT"


class Dtype(StrEnum):
    # Types of dose-response datasets
    DICHOTOMOUS = "D"
    CONTINUOUS = "C"
    CONTINUOUS_INDIVIDUAL = "CI"
    NESTED_DICHOTOMOUS = "ND"

    @classmethod
    def CONTINUOUS_DTYPES(cls):
        return {cls.CONTINUOUS, cls.CONTINUOUS_INDIVIDUAL}


# model names
class Models(StrEnum):
    Weibull = "Weibull"
    LogProbit = "LogProbit"
    Probit = "Probit"
    QuantalLinear = "Quantal Linear"
    Multistage = "Multistage"
    Gamma = "Gamma"
    Logistic = "Logistic"
    LogLogistic = "LogLogistic"
    DichotomousHill = "Dichotomous-Hill"
    Linear = "Linear"
    Polynomial = "Polynomial"
    Power = "Power"
    Exponential = "Exponential"
    ExponentialM2 = "Exponential-M2"
    ExponentialM3 = "Exponential-M3"
    ExponentialM4 = "Exponential-M4"
    ExponentialM5 = "Exponential-M5"
    Hill = "Hill"
    MultiplicativeHill = "Multiplicative Hill"
    InverseExponential = "Inverse Exponential"
    Lognormal = "Lognormal"
    ContinuousGamma = "Continuous Gamma"
    LMS2 = "LMS 2-Stage"
    NestedLogistic = "Nested Logistic"
    NCTR = "NCTR"

    @classmethod
    def VARIABLE_POLYNOMIAL(cls):
        return {cls.Multistage, cls.Polynomial}


class LogicBin(IntEnum):
    NO_CHANGE = 0, "valid", "✓", "Viable"
    WARNING = 1, "warning", "?", "Questionable"
    FAILURE = 2, "failure", "✕", "Unusable"

    def __new__(cls, value, text, icon, bin):
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.text = text
        obj.icon = icon
        obj.bin = bin
        return obj


class BmdModelSchema(BaseModel):
    id: int
    verbose: str
    bmds_model_form_str: str = Field(alias="model_form_str")


class DichotomousModel(BmdModelSchema):
    params: tuple[str, ...]

    @property
    def num_params(self):
        return len(self.params)


class DichotomousModelChoices(Enum):
    hill = DichotomousModel(
        id=dich_model.d_hill.value,
        verbose="Hill",
        params=("g", "v", "a", "b"),
        model_form_str="P[dose] = g + (v - v * g) / (1 + exp(-a - b * Log(dose)))",
    )
    gamma = DichotomousModel(
        id=dich_model.d_gamma.value,
        verbose="Gamma",
        params=("g", "a", "b"),
        model_form_str="P[dose]= g + (1 - g) * CumGamma(b * dose, a)",
    )
    logistic = DichotomousModel(
        id=dich_model.d_logistic.value,
        verbose="Logistic",
        params=("a", "b"),
        model_form_str="P[dose] = 1 / [1 + exp(-a - b * dose)]",
    )
    loglogistic = DichotomousModel(
        id=dich_model.d_loglogistic.value,
        verbose="LogLogistic",
        params=("g", "a", "b"),
        model_form_str="P[dose] = g + (1 - g)/(1 + exp(-a - b * Log(dose)))",
    )
    logprobit = DichotomousModel(
        id=dich_model.d_logprobit.value,
        verbose="LogProbit",
        params=("g", "a", "b"),
        model_form_str="P[dose] = g + (1 - g) * CumNorm(a + b * Log(Dose))",
    )
    multistage = DichotomousModel(
        id=dich_model.d_multistage.value,
        verbose="Multistage",
        params=("g", "b1", "bN"),
        model_form_str="P[dose] = g + (1 - g) * (1 - exp(-b1 * dose^1 - b2 * dose^2 - ...))",
    )
    probit = DichotomousModel(
        id=dich_model.d_probit.value,
        verbose="Probit",
        params=("a", "b"),
        model_form_str="P[dose] = CumNorm(a + b * Dose)",
    )
    qlinear = DichotomousModel(
        id=dich_model.d_qlinear.value,
        verbose="Quantal Linear",
        params=("g", "b"),
        model_form_str="P[dose] = g + (1 - g) * (1 - exp(-b * dose)",
    )
    weibull = DichotomousModel(
        id=dich_model.d_weibull.value,
        verbose="Weibull",
        params=("g", "a", "b"),
        model_form_str="P[dose] = g + (1 - g) * (1 - exp(-b * dose^a))",
    )


class ContinuousModel(BmdModelSchema):
    params: tuple[str, ...]
    variance_params: tuple[str, ...]


class ContinuousModelChoices(Enum):
    power = ContinuousModel(
        id=cont_model.power.value,
        verbose="Power",
        params=("g", "v", "n"),
        variance_params=("rho", "alpha"),
        model_form_str="P[dose] = g + v * dose ^ n",
    )
    hill = ContinuousModel(
        id=cont_model.hill.value,
        verbose="Hill",
        params=("g", "v", "k", "n"),
        variance_params=("rho", "alpha"),
        model_form_str="P[dose] = g + v * dose ^ n / (k ^ n + dose ^ n)",
    )
    polynomial = ContinuousModel(
        id=cont_model.polynomial.value,
        verbose="Polynomial",
        params=("g", "b1", "bN"),
        variance_params=("rho", "alpha"),
        model_form_str="P[dose] = g + b1*dose + b2*dose^2 + b3*dose^3...",
    )
    exp_m3 = ContinuousModel(
        id=cont_model.exp_3.value,
        verbose="Exponential 3",
        params=("a", "b", "c", "d"),
        variance_params=("rho", "log-alpha"),
        model_form_str="P[dose] = a * exp(±1 * (b * dose) ^ d)",
    )
    exp_m5 = ContinuousModel(
        id=cont_model.exp_5.value,
        verbose="Exponential 5",
        params=("a", "b", "c", "d"),
        variance_params=("rho", "log-alpha"),
        model_form_str="P[dose] = a * (c - (c - 1) * exp(-(b * dose) ^ d)",
    )
    mult_hill = ContinuousModel(
        id=loud_cont_model.l_hill_efsa.value,
        verbose="Multiplicative Hill",
        params=("a", "b", "c", "d"),
        variance_params=("rho", "log-alpha"),
        model_form_str="P[dose] = a * (1 + (c - 1) * dose ^ d / (b ^ d + dose ^ d))",
        )
    inverse_exp = ContinuousModel(
        id=loud_cont_model.l_invexp_efsa.value,
        verbose="Inverse Exponential",
        params=("a", "b", "c", "d"),
        variance_params=("rho", "log-alpha"),
        model_form_str="P[dose] = a * (1 + (c - 1) * exp(-b * dose ^ -d))",
        )
    lognormal = ContinuousModel(
        id=loud_cont_model.l_lognormal_efsa.value,
        verbose="Lognormal",
        params=("a", "b", "c", "d"),
        variance_params=("rho", "log-alpha"),
        model_form_str="P[dose] = a * (1 + (c - 1) * CumNorm(Log(b) + d * Log(dose)))",
        )
    cont_gamma = ContinuousModel(
        id=loud_cont_model.l_gamma_efsa.value,
        verbose="Continuous Gamma",
        params=("a", "b", "c", "d"),
        variance_params=("rho", "log-alpha"),
        model_form_str="P[dose] = a * (1 + (c - 1) * CumGamma(b * dose, d))",
        )
    lms = ContinuousModel(
        id=loud_cont_model.l_lms_efsa.value,
        verbose="LMS 2-Stage",
        params=("a", "b", "c", "d"),
        variance_params=("rho", "log-alpha"),
        model_form_str="P[dose] = a * (1 + (c - 1) * (1 - exp(-b * dose - d * dose ^ 2)))",
        )


class DistType(IntEnum):
    normal = 1, "Normal", "Constant variance"  # f(i) = a * x(i)
    normal_ncv = 2, "Normal", "Nonconstant variance"  # f(i) = a * x(i) ^ p
    log_normal = 3, "Lognormal", "Constant variance"

    def __new__(cls, value, distribution_type, variance_model):
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.distribution_type = distribution_type
        obj.variance_model = variance_model
        return obj


class PriorDistribution(IntEnum):
    Uniform = 0
    Normal = 1
    Lognormal = 2


class PriorClass(IntEnum):
    frequentist_unrestricted = 0, "Frequentist unrestricted", "Unrestricted", False
    frequentist_restricted = 1, "Frequentist restricted", "Restricted", False
    bayesian = 2, "Bayesian", "N/A", True

    def __new__(cls, value, label, restriction, is_bayesian):
        obj = int.__new__(cls, value)
        obj._value_ = value
        obj.label = label
        obj.restriction = restriction
        obj.is_bayesian = is_bayesian
        return obj


class NestedDichotomousModel(BmdModelSchema):
    params: tuple[str, ...]


class NestedDichotomousModelChoices(Enum):
    logistic = NestedDichotomousModel(
        id=nested_model.nlogistic.value,
        verbose="Nested Logistic",
        params=("a", "b", "theta1", "theta2", "rho", "phi"),
        model_form_str="P[dose] = a + θ₁ * r + (1 - a - θ₁ * r) / (1 + e^[-b - θ₂ * r - p * ln(dose)])",
    )
    nctr = NestedDichotomousModel(
        id=nested_model.nctr.value,
        verbose="NCTR",
        params=("a", "b", "theta1", "theta2", "rho", "phi"),
        model_form_str="P[dose] =  1 - e^[-(a + θ₁ * (r - rm)) - (b + θ₂ * (r - rm)) * dose^p]",
    )
