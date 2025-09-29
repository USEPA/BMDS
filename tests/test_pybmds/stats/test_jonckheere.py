import numpy as np
import pytest

from pybmds.stats.jonckheere import jonckheere


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

        # non-unique x
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
