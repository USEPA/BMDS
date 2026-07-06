import numpy as np
import pytest

from pybmds.constants import Dtype
from pybmds.models.base import BmdModelAveragingSchema, BmdModelSchema, cdf_df


def test_cdf_df_combines_chains_for_empirical_cdf():
    cdf = cdf_df(np.array([[4.0, 1.0, np.nan], [3.0, 2.0, np.inf]]), n_points=3)

    assert cdf["Percentile"].to_list() == pytest.approx([0.5, 50.0, 99.5])
    assert cdf["BMD"].to_list() == pytest.approx(
        np.percentile([1.0, 2.0, 3.0, 4.0], [0.5, 50, 99.5])
    )


def test_cdf_df_preserves_bmds_cdf_format():
    cdf = cdf_df(np.array([[10.0, 20.0, 30.0], [0.1, 0.5, 0.9]]))

    assert cdf["Percentile"].to_list() == pytest.approx([10.0, 50.0, 90.0])
    assert cdf["BMD"].to_list() == pytest.approx([10.0, 20.0, 30.0])


def test_cdf_df_handles_empty_draws_and_rejects_higher_dimensions():
    empty = cdf_df(np.array([np.nan, np.inf]))
    assert empty.empty
    assert empty.columns.to_list() == ["Percentile", "BMD"]

    with pytest.raises(ValueError, match="Unsupported CDF input shape"):
        cdf_df(np.zeros((2, 2, 2)))


def test_model_schema_subclass_dispatch_and_invalid_dtypes():
    assert BmdModelSchema.get_subclass(Dtype.DICHOTOMOUS).__name__.endswith("DichotomousSchema")
    assert BmdModelSchema.get_subclass(Dtype.CONTINUOUS).__name__.endswith("ContinuousSchema")
    assert BmdModelSchema.get_subclass(Dtype.NESTED_DICHOTOMOUS).__name__.endswith(
        "NestedDichotomousSchema"
    )
    assert BmdModelAveragingSchema.get_subclass(Dtype.DICHOTOMOUS).__name__.endswith(
        "DichotomousSchema"
    )
    assert BmdModelAveragingSchema.get_subclass(Dtype.CONTINUOUS).__name__.endswith(
        "ContinuousSchema"
    )
    with pytest.raises(ValueError, match="Invalid dtype"):
        BmdModelSchema.get_subclass("invalid")
    with pytest.raises(ValueError, match="Invalid dtype"):
        BmdModelAveragingSchema.get_subclass(Dtype.NESTED_DICHOTOMOUS)
