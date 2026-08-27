import numpy as np
import pytest

from pybmds.stats.jonckheere import Approach, jonckheere


@pytest.fixture
def valid_x():
    return np.linspace(1, 8, 8)


@pytest.fixture
def rep_x():
    return np.array([1.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])


@pytest.fixture
def valid_group():
    return np.repeat(np.linspace(1, 4, 4), 2)


class TestJonckheere:
    def test_obs_warning(self, rep_x, valid_group):
        with pytest.warns(UserWarning, match="total observations < 30"):
            jonckheere(rep_x, valid_group, hypothesis="increasing")

    def test_exact_result(self, valid_x, valid_group):
        result = jonckheere(valid_x, valid_group, hypothesis="increasing")
        assert pytest.approx(result.statistic) == 24.0
        assert pytest.approx(result.p_value, abs=1e-3) == 0.00109

        result = jonckheere(valid_x, valid_group, hypothesis="two-sided")
        assert pytest.approx(result.statistic) == 24.0
        assert pytest.approx(result.p_value, abs=1e-4) == 0.000793

    def test_valid_result(self, valid_x, valid_group):
        result = jonckheere(valid_x, valid_group, hypothesis="increasing")
        assert 0 <= result.p_value <= 1

    def test_result_table(self, valid_x, valid_group):
        result = jonckheere(valid_x, valid_group)
        assert isinstance(result.tbl(), str)

    @pytest.mark.parametrize("hypothesis", ("increasing", "decreasing", "two-sided"))
    def test_hypothesis_paths(self, hypothesis, valid_x, valid_group):
        # run all code paths and confirm we always get a p_value.

        # unique x
        result = jonckheere(valid_x, valid_group, hypothesis=hypothesis)
        assert 0 <= result.p_value <= 1

        with pytest.warns(
            UserWarning,
            match="P-value estimated using normal distribution; total observations < 30",
        ):
            result = jonckheere(
                np.repeat(np.linspace(1, 10, 10), 2),
                np.repeat(np.linspace(1, 5, 5), 4),
                hypothesis=hypothesis,
            )

        assert 0 <= result.p_value <= 1

        # permutations
        result = jonckheere(valid_x, valid_group, hypothesis=hypothesis, nperm=10)
        assert 0 <= result.p_value <= 1

    def test_data_not_numeric(self, valid_group):
        x = np.array("a b c".split())
        with pytest.raises(ValueError, match="Data needs to be numeric"):
            jonckheere(x, valid_group, hypothesis="two-sided")

    def test_group_not_numeric(self, valid_x):
        group = np.array("a b c".split())
        with pytest.raises(ValueError, match="Group needs to be numeric or ordered factor"):
            jonckheere(valid_x, group, hypothesis="two-sided")

    def test_group_data_different_lengths(self, valid_x):
        group = np.array([1, 1, 2, 2, 3])
        with pytest.raises(ValueError, match="Data and group values need to be the same length"):
            jonckheere(valid_x, group, hypothesis="two-sided")

    def test_hypothesis_not_valid(self, valid_x, valid_group):
        hypothesis = "test"
        with pytest.raises(ValueError, match="'test' is not a valid Hypothesis"):
            jonckheere(valid_x, valid_group, hypothesis)

    def test_data_empty(self, valid_group):
        x = np.repeat(np.nan, 8)
        with pytest.raises(
            ValueError, match="Either data or group is missing for all observations"
        ):
            jonckheere(x, valid_group, hypothesis="two-sided")

    def test_group_empty(self, valid_x):
        group = np.repeat(np.nan, 8)
        with pytest.raises(
            ValueError, match="Either data or group is missing for all observations"
        ):
            jonckheere(valid_x, group, hypothesis="two-sided")

    def test_one_group(self, valid_x):
        group = np.ones(8)
        with pytest.raises(ValueError, match="Only one group has non-missing data"):
            jonckheere(valid_x, group, hypothesis="two-sided")

    def test_exact_increasing_trend_small_pvalue(self):
        """Exact method should give a small p-value for a strong increasing trend."""
        x = np.array(
            [
                1,
                2,
                3,
                10,
                11,
                12,
                20,
                21,
                22,
            ],
            dtype=float,
        )

        group = np.array(
            [
                0,
                0,
                0,
                1,
                1,
                1,
                2,
                2,
                2,
            ],
            dtype=float,
        )

        result = jonckheere(x, group, hypothesis="increasing")

        assert result.approach == Approach.exact
        assert result.p_value < 0.05

    def test_exact_decreasing_trend_small_pvalue(self):
        """Exact method should give a small p-value for a strong decreasing trend."""
        x = np.array(
            [
                20,
                21,
                22,
                10,
                11,
                12,
                1,
                2,
                3,
            ],
            dtype=float,
        )

        group = np.array(
            [
                0,
                0,
                0,
                1,
                1,
                1,
                2,
                2,
                2,
            ],
            dtype=float,
        )

        result = jonckheere(x, group, hypothesis="decreasing")

        assert result.approach == Approach.exact
        assert result.p_value < 0.05

    def test_permutation_increasing_trend_small_pvalue(self):
        """Permutation method should use the upper tail for an increasing trend."""
        x = np.array(
            [
                1,
                2,
                3,
                10,
                11,
                12,
                20,
                21,
                22,
            ],
            dtype=float,
        )

        group = np.array(
            [
                0,
                0,
                0,
                1,
                1,
                1,
                2,
                2,
                2,
            ],
            dtype=float,
        )

        result = jonckheere(
            x,
            group,
            hypothesis="increasing",
            nperm=10000,
            seed=123,
        )

        assert result.approach == Approach.permutation
        assert result.p_value < 0.05

    def test_permutation_decreasing_trend_small_pvalue(self):
        """Permutation method should use the lower tail for a decreasing trend."""
        x = np.array(
            [
                20,
                21,
                22,
                10,
                11,
                12,
                1,
                2,
                3,
            ],
            dtype=float,
        )

        group = np.array(
            [
                0,
                0,
                0,
                1,
                1,
                1,
                2,
                2,
                2,
            ],
            dtype=float,
        )

        result = jonckheere(
            x,
            group,
            hypothesis="decreasing",
            nperm=10000,
            seed=123,
        )

        assert result.approach == Approach.permutation
        assert result.p_value < 0.05
