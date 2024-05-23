import math
from typing import Any, ClassVar, Self

from pydantic import BaseModel

from .. import constants
from ..constants import BMDS_BLANK_VALUE, DistType
from .constants import RuleClass

Number = float | int


class CheckResponse(BaseModel):
    logic_bin: constants.LogicBin
    message: str

    @classmethod
    def success(cls) -> Self:
        return CheckResponse(logic_bin=constants.LogicBin.NO_CHANGE, message="")


class Check:
    _enabled_attribute: ClassVar[dict] = {
        constants.Dtype.DICHOTOMOUS: "enabled_dichotomous",
        constants.Dtype.CONTINUOUS: "enabled_continuous",
        constants.Dtype.CONTINUOUS_INDIVIDUAL: "enabled_continuous",
        constants.Dtype.NESTED_DICHOTOMOUS: "enabled_nested",
    }

    @classmethod
    def check(cls, dataset, model, rule_settings) -> CheckResponse:
        is_enabled = getattr(rule_settings, cls._enabled_attribute[dataset.dtype])
        if not is_enabled:
            return CheckResponse.success()
        if failure_msg := cls.run_check(dataset, model, rule_settings):
            return CheckResponse(
                logic_bin=rule_settings.failure_bin,
                message=failure_msg,
            )
        return CheckResponse.success()

    @classmethod
    def run_check(cls, dataset, model, rule_settings) -> str | None:
        """Execute a check; return an error message if failure, else None"""
        raise NotImplementedError("")


def is_valid_number(value: Any) -> bool:
    try:
        return math.isfinite(value) and value != BMDS_BLANK_VALUE
    except TypeError:
        return False


def number_or_none(value: Any) -> Number | None:
    if not is_valid_number(value):
        return None
    return value


def is_gte(value: Any, threshold: float) -> bool:
    if not is_valid_number(value):
        return True  # if value is invalid; other checks will report error
    return value >= threshold


def is_lte(value: Any, threshold: float) -> bool:
    if not is_valid_number(value):
        return True  # if value is invalid; other checks will report error
    return value <= threshold


# existence checks
# --------------------------------------------------------------------------------------------------
class ExistenceCheck(Check):
    failure_message_name: str

    @classmethod
    def get_value(cls, dataset, model) -> Number | None: ...

    @classmethod
    def run_check(cls, dataset, model, rule_settings) -> str | None:
        value = cls.get_value(dataset, model)
        if not is_valid_number(value):
            return f"{cls.failure_message_name} does not exist"


class AicExists(ExistenceCheck):
    failure_message_name = "AIC"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        return model.results.get_parameter("aic")


class BmdExists(ExistenceCheck):
    failure_message_name = "BMD"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        return model.results.get_parameter("bmd")


class RoiExists(ExistenceCheck):
    failure_message_name = "Residual of interest"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        return model.results.get_parameter("roi")


class BmdlExists(ExistenceCheck):
    failure_message_name = "BMDL"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        # bmdl must also be non-zero
        bmdl = model.results.get_parameter("bmdl")
        if bmdl == 0:
            return BMDS_BLANK_VALUE
        return bmdl


class BmduExists(ExistenceCheck):
    failure_message_name = "BMDU"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        return model.results.get_parameter("bmdu")


# greater than threshold checks
# --------------------------------------------------------------------------------------------------
class ShouldBeGreaterThan(Check):
    failure_message_name: str

    @classmethod
    def get_value(cls, dataset, model) -> Number | None: ...

    @classmethod
    def run_check(cls, dataset, model, rule_settings) -> str | None:
        value = cls.get_value(dataset, model)
        threshold = rule_settings.threshold
        if not is_gte(value, threshold):
            return f"{cls.failure_message_name} < {threshold}"


class GoodnessOfFit(ShouldBeGreaterThan):
    failure_message_name = "Goodness of fit p-value"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        # only run test if DOF is > 0
        dof = model.results.get_parameter("dof")
        if dof <= constants.ZEROISH:
            return None
        return model.results.get_parameter("pvalue")


class GoodnessOfFitCancer(GoodnessOfFit):
    failure_message_name = "Goodness of fit p-value (cancer)"


# less than threshold checks
# --------------------------------------------------------------------------------------------------
class ShouldBeLessThan(Check):
    failure_message_name: str

    @classmethod
    def get_value(cls, dataset, model) -> Number | None: ...

    @classmethod
    def run_check(cls, dataset, model, rule_settings) -> str | None:
        value = cls.get_value(dataset, model)
        threshold = rule_settings.threshold
        if not is_lte(value, threshold):
            return f"{cls.failure_message_name} > {threshold}"


class LargeRoi(ShouldBeLessThan):
    failure_message_name = "|Residual near BMD|"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        roi = model.results.get_parameter("roi")
        if is_valid_number(roi):
            return abs(roi)
        return None


class BmdBmdlRatio(ShouldBeLessThan):
    failure_message_name = "BMD/BMDL ratio"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        bmd = model.results.get_parameter("bmd")
        bmdl = model.results.get_parameter("bmdl")
        if is_valid_number(bmd) and is_valid_number(bmdl) and bmdl > 0:
            return bmd / bmdl
        return None


class LowBmd(ShouldBeLessThan):
    failure_message_name = "lowest dose/BMD ratio"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        min_dose = min([dose for dose in dataset.doses if dose > 0])
        bmd = model.results.get_parameter("bmd")
        if is_valid_number(min_dose) and is_valid_number(bmd) and bmd > 0:
            return min_dose / float(bmd)
        return None


class LowBmdl(ShouldBeLessThan):
    failure_message_name = "lowest dose/BMDL ratio"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        min_dose = min([d for d in dataset.doses if d > 0])
        bmdl = model.results.get_parameter("bmdl")
        if is_valid_number(min_dose) and is_valid_number(bmdl) and bmdl > 0:
            return min_dose / float(bmdl)
        return None


class HighBmd(ShouldBeLessThan):
    failure_message_name = "BMD/highest dose ratio"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        max_dose = max(dataset.doses)
        bmd = model.results.get_parameter("bmd")
        if is_valid_number(max_dose) and is_valid_number(bmd) and bmd > 0:
            return bmd / float(max_dose)
        return None


class HighBmdl(ShouldBeLessThan):
    failure_message_name = "BMDL/highest dose ratio"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        max_dose = max(dataset.doses)
        bmdl = model.results.get_parameter("bmdl")
        if is_valid_number(max_dose) and is_valid_number(bmdl) and bmdl > 0:
            return bmdl / float(max_dose)
        return None


class HighControlResidual(ShouldBeLessThan):
    failure_message_name = "Residual at control"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        return abs(model.results.get_parameter("roi_control"))


class ControlStdevFit(ShouldBeLessThan):
    failure_message_name = "Control stdev. fit"

    @classmethod
    def get_value(cls, dataset, model) -> Number | None:
        modeled = model.results.gof.est_sd[0]
        actual = model.results.gof.calc_sd[0]
        return modeled / max(actual, 1e-6)


# assorted checks
# --------------------------------------------------------------------------------------------------
class VarianceFit(Check):
    @classmethod
    def run_check(cls, dataset, model, rule_settings) -> str | None:
        constant_variance = model.settings.disttype != DistType.normal_ncv
        test = 2 if constant_variance else 3
        model_str = "Constant" if constant_variance else "Nonconstant"
        pvalue = model.results.tests.p_values[test - 1]
        if is_valid_number(pvalue) and pvalue < rule_settings.threshold:
            return f"{model_str} variance test failed (Test {test} p-value < {rule_settings.threshold})"


class VarianceType(Check):
    @classmethod
    def run_check(cls, dataset, model, rule_settings) -> str | None:
        constant_variance = model.settings.disttype != DistType.normal_ncv
        p_value2 = model.results.tests.p_values[1]
        if is_valid_number(p_value2):
            if constant_variance and p_value2 < rule_settings.threshold:
                return f"Incorrect variance model (Test 2 p-value < {rule_settings.threshold})"
            if not constant_variance and p_value2 > rule_settings.threshold:
                return f"Incorrect variance model (Test 2 p-value > {rule_settings.threshold})"


class NoDegreesOfFreedom(Check):
    @classmethod
    def run_check(cls, dataset, model, rule_settings) -> str | None:
        value = model.results.get_parameter("dof")
        if value <= constants.ZEROISH:
            return "Zero degrees of freedom; saturated model"


class Warnings(Check):
    @classmethod
    def run_check(cls, settings, dataset, output) -> str | None:
        return None


RULE_MAP = {
    RuleClass.gof: GoodnessOfFit,
    RuleClass.dof_zero: NoDegreesOfFreedom,
    RuleClass.high_bmd: HighBmd,
    RuleClass.warnings: Warnings,
    RuleClass.high_bmdl: HighBmdl,
    RuleClass.roi_large: LargeRoi,
    RuleClass.gof_cancer: GoodnessOfFitCancer,
    RuleClass.aic_missing: AicExists,
    RuleClass.bmd_missing: BmdExists,
    RuleClass.roi_missing: RoiExists,
    RuleClass.bmdl_missing: BmdlExists,
    RuleClass.bmdu_missing: BmduExists,
    RuleClass.low_bmd_fail: LowBmd,
    RuleClass.low_bmd_warn: LowBmd,
    RuleClass.variance_fit: VarianceFit,
    RuleClass.low_bmdl_fail: LowBmdl,
    RuleClass.low_bmdl_warn: LowBmdl,
    RuleClass.variance_type: VarianceType,
    RuleClass.control_stdev_fit: ControlStdevFit,
    RuleClass.bmd_bmdl_ratio_fail: BmdBmdlRatio,
    RuleClass.bmd_bmdl_ratio_warn: BmdBmdlRatio,
    RuleClass.control_residual_high: HighControlResidual,
}
