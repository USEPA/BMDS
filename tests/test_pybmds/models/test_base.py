import numpy as np
import pytest

import pybmds
from pybmds.constants import BMDS_BLANK_VALUE, Dtype
from pybmds.models.base import BmdModelAveragingSchema, BmdModelSchema, cdf_df


def test_cdf_df_combines_chains_for_empirical_cdf():
    cdf = cdf_df(
        np.array([[4.0, 1.0, np.nan], [3.0, 2.0, np.inf], [BMDS_BLANK_VALUE, 5.0, 6.0]]),
        n_points=3,
    )

    assert cdf["Percentile"].to_list() == pytest.approx([0.5, 50.0, 99.5])
    assert cdf["BMD"].to_list() == pytest.approx(
        np.percentile([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [0.5, 50, 99.5])
    )


def test_cdf_df_preserves_bmds_cdf_format():
    cdf = cdf_df(np.array([[10.0, BMDS_BLANK_VALUE, 30.0], [0.1, 0.5, 0.9]]))

    assert cdf["Percentile"].to_list() == pytest.approx([10.0, 90.0])
    assert cdf["BMD"].to_list() == pytest.approx([10.0, 30.0])


def test_cdf_df_handles_empty_draws_and_rejects_higher_dimensions():
    empty = cdf_df(np.array([np.nan, np.inf, BMDS_BLANK_VALUE]))
    assert empty.empty
    assert empty.columns.to_list() == ["Percentile", "BMD"]

    empty_from_cdf_shape = cdf_df(np.array([[BMDS_BLANK_VALUE, np.nan], [0.25, 0.75]]))
    assert empty_from_cdf_shape["BMD"].to_list() == pytest.approx(
        np.percentile([0.25, 0.75], empty_from_cdf_shape["Percentile"])
    )

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


def test_plot_handles_completed_model_without_bmd_lines(cdataset2):
    model = pybmds.models.continuous.Power(cdataset2)
    model.execute()
    model.results = model.results.model_copy(
        update={
            "bmd": BMDS_BLANK_VALUE,
            "bmdl": BMDS_BLANK_VALUE,
            "bmdu": BMDS_BLANK_VALUE,
        }
    )

    fig = model.plot()

    labels = [text.get_text() for text in fig.gca().get_legend().get_texts()]
    assert "BMDL-BMD-BMDU" not in labels
