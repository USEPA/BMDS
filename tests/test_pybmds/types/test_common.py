import numpy as np
import pytest
from pydantic import BaseModel, ValidationError

from pybmds.types.common import NumpyFloatArray, residual_of_interest


def test_residual_of_interest():
    # failure case
    assert residual_of_interest(bmd=-9999, doses=[1, 2, 3], residuals=[0.1, 0.2, 0.3]) == -9999

    # above value
    assert residual_of_interest(bmd=2.1, doses=[1, 2, 3], residuals=[0.1, 0.2, 0.3]) == 0.2

    # below value
    assert residual_of_interest(bmd=1.9, doses=[1, 2, 3], residuals=[0.1, 0.2, 0.3]) == 0.2

    # exactly between value; chooses smaller value
    assert residual_of_interest(bmd=1, doses=[0, 2, 4], residuals=[0.1, 0.2, 0.3]) == 0.1


class ExampleModel(BaseModel):
    d: NumpyFloatArray


class TestNumpyFloatArray:
    def test_successes(self):
        for data in [[1], [1, 2, 3], [[1, 2], [3, 4]]]:
            model = ExampleModel(d=data)
            assert np.allclose(model.d, data)
            assert model.d.dtype == float
            assert isinstance(model.d, np.ndarray)

    def test_failures(self):
        for data in ["a", ["a", "a"], None]:
            with pytest.raises(ValidationError):
                ExampleModel(d=data)
