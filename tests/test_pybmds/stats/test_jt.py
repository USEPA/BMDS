import numpy as np
import pytest

from pybmds.stats.jt import JonckheereResult, jonckheere

dummy_data = np.linspace(1, 8, 8)
dummy_group = np.repeat(np.linspace(1, 4, 4), 2)
dummy_uneven = np.array([1, 1, 2, 2, 3])
dummy_nonnumeric = np.array(["a", "b", "c", "d", "e", "f", "g", "h"])
dummy_ones = np.ones(8)
dummy_empty = np.repeat(np.nan, 8)


class TestJonckheereTrendTest:
    def test_valid_result(self):
        x = dummy_data
        group = dummy_group
        result = jonckheere(x, group, alternative="increasing")
        assert isinstance(result, JonckheereResult)
        assert 0 <= result.p_value <= 1

    def test_data_not_numeric(self):
        x = dummy_nonnumeric
        group = dummy_group
        with pytest.raises(ValueError, match="Data needs to be numeric"):
            jonckheere(x, group, alternative="two-sided")

    def test_group_not_numeric(self):
        x = dummy_data
        group = dummy_nonnumeric
        with pytest.raises(ValueError, match="Group needs to be numeric or ordered factor"):
            jonckheere(x, group, alternative="two-sided")

    def test_group_data_different_lengths(self):
        x = dummy_data
        group = dummy_uneven
        with pytest.raises(ValueError, match="Data and group values need to be the same length"):
            jonckheere(x, group, alternative="two-sided")

    def test_alternative_not_valid(self):
        x = dummy_data
        group = dummy_group
        alternative = "test"
        with pytest.raises(ValueError, match="'test' is not a valid Alternative"):
            jonckheere(x, group, alternative)

    def test_data_empty(self):
        x = dummy_empty
        group = dummy_group
        with pytest.raises(
            ValueError, match="Either data or group is missing for all observations"
        ):
            jonckheere(x, group, alternative="two-sided")

    def test_group_empty(self):
        x = dummy_data
        group = dummy_empty
        with pytest.raises(
            ValueError, match="Either data or group is missing for all observations"
        ):
            jonckheere(x, group, alternative="two-sided")

    def test_one_group(self):
        x = dummy_data
        group = dummy_ones
        with pytest.raises(ValueError, match="Only one group has non-missing data"):
            jonckheere(x, group, alternative="two-sided")
