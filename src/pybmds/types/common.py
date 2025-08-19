from typing import Any

import numpy as np
from pydantic import GetCoreSchemaHandler
from pydantic_core import CoreSchema, core_schema

from .. import bmdscore
from ..constants import BMDS_BLANK_VALUE


def residual_of_interest(bmd: float, doses: list[float], residuals: list[float]) -> float:
    if bmd <= 0:
        return BMDS_BLANK_VALUE
    diffs = [abs(bmd - dose) for dose in doses]
    index = diffs.index(min(diffs))
    return residuals[index]


def clean_array(arr: np.ndarray) -> np.ndarray:
    return np.nan_to_num(
        arr, nan=BMDS_BLANK_VALUE, posinf=BMDS_BLANK_VALUE, neginf=BMDS_BLANK_VALUE
    )


class PydanticNumpyArray(np.ndarray):
    # pydantic friendly numpy arrays

    @classmethod
    def __get_pydantic_core_schema__(
        cls, source_type: Any, handler: GetCoreSchemaHandler
    ) -> CoreSchema:
        return core_schema.no_info_after_validator_function(
            cls.validate,
            handler("list"),
            serialization=core_schema.plain_serializer_function_ser_schema(lambda x: x.tolist()),
        )


class NumpyIntArray(PydanticNumpyArray):
    @classmethod
    def validate(cls, v):
        try:
            return np.asarray(v, dtype="int")
        except TypeError as err:
            raise ValueError("invalid np.ndarray format") from err


class NumpyFloatArray(PydanticNumpyArray):
    # Numpy arrays, agumented
    @classmethod
    def validate(cls, v):
        try:
            return np.asarray(v, dtype="float")
        except TypeError as err:
            raise ValueError("invalid np.ndarray format") from err


def inspect_cpp_obj(obj: Any) -> str:
    """Prints struct contents using C++ function

    Captures stream from C++ function and stores in list of strings.

    Args:
        obj (Any): the object to inspect
    """
    return bmdscore.print(obj, False)


BOUND_FOOTNOTE = """
Standard errors estimates are not generated for parameters estimated on corresponding bounds,
although sampling error is present for all parameters, as a rule. Standard error estimates may not
be reliable as a basis for confidence intervals or tests when one or more parameters are on bounds.
"""


CONTINUOUS_TEST_FOOTNOTES = """\
Test 1: Test the null hypothesis that responses and variances don't differ among dose levels
(A2 vs R).  If this test fails to reject the null hypothesis (p-value > 0.05), there may not be
a dose-response.

Test 2: Test the null hypothesis that variances are homogenous (A1 vs A2).  If this test fails to
reject the null hypothesis (p-value > 0.05), the simpler constant variance model may be appropriate.

Test 3: Test the null hypothesis that the variances are adequately modeled (A3 vs A2). If this test
fails to reject the null hypothesis (p-value > 0.05), it may be inferred that the variances have
been modeled appropriately.

Test 4: Test the null hypothesis that the model for the mean fits the data (Fitted vs A3). If this
test fails to reject the null hypothesis (p-value > 0.1), the user has support for use of the
selected model.
"""
