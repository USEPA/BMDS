import numpy as np
import pytest

from pybmds.stats.jt import JonckheereResult, jonckheere


@pytest.fixture
def valid_x():
    return np.linspace(1, 8, 8)


@pytest.fixture
def valid_group():
    return np.repeat(np.linspace(1, 4, 4), 2)


class TestJonckheere:
    def test_valid_result(self, valid_x, valid_group):
        result = jonckheere(valid_x, valid_group, alternative="increasing")
        assert isinstance(result, JonckheereResult)
        assert 0 <= result.p_value <= 1

    @pytest.mark.parametrize("alternative", ("increasing", "decreasing", "two-sided"))
    def test_alternative_paths(self, alternative, valid_x, valid_group):
        # run all code paths and confirm we always get a p_value.

        # unique x
        result = jonckheere(valid_x, valid_group, alternative=alternative)
        assert 0 <= result.p_value <= 1

        # non-unique x
        result = jonckheere(
            np.repeat(np.linspace(1, 10, 10), 2),
            np.repeat(np.linspace(1, 5, 5), 4),
            alternative=alternative,
        )
        assert 0 <= result.p_value <= 1

        # permutations
        result = jonckheere(valid_x, valid_group, alternative=alternative, nperm=10)
        assert 0 <= result.p_value <= 1

    def test_data_not_numeric(self, valid_group):
        x = np.array("a b c".split())
        with pytest.raises(ValueError, match="Data needs to be numeric"):
            jonckheere(x, valid_group, alternative="two-sided")

    def test_group_not_numeric(self, valid_x):
        group = np.array("a b c".split())
        with pytest.raises(ValueError, match="Group needs to be numeric or ordered factor"):
            jonckheere(valid_x, group, alternative="two-sided")

    def test_group_data_different_lengths(self, valid_x):
        group = np.array([1, 1, 2, 2, 3])
        with pytest.raises(ValueError, match="Data and group values need to be the same length"):
            jonckheere(valid_x, group, alternative="two-sided")

    def test_alternative_not_valid(self, valid_x, valid_group):
        alternative = "test"
        with pytest.raises(ValueError, match="'test' is not a valid Alternative"):
            jonckheere(valid_x, valid_group, alternative)

    def test_data_empty(self, valid_group):
        x = np.repeat(np.nan, 8)
        with pytest.raises(
            ValueError, match="Either data or group is missing for all observations"
        ):
            jonckheere(x, valid_group, alternative="two-sided")

    def test_group_empty(self, valid_x):
        group = np.repeat(np.nan, 8)
        with pytest.raises(
            ValueError, match="Either data or group is missing for all observations"
        ):
            jonckheere(valid_x, group, alternative="two-sided")

    def test_one_group(self, valid_x):
        group = np.ones(8)
        with pytest.raises(ValueError, match="Only one group has non-missing data"):
            jonckheere(valid_x, group, alternative="two-sided")
