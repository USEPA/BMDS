import numpy as np
import pytest

from pybmds.models.base import cdf_df


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
