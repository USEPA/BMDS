import numpy as np
import pytest

from pybmds.stats.jt import jonckheere

dummy_data = np.array([1, 2, 3, 4, 5, 6, 7, 8])
dummy_group = np.array([1, 1, 2, 2, 3, 3, 4, 4])
dummy_uneven = np.array([1, 1, 2, 2, 3])
alternative = "two.sided"
dummy_nonnumeric = np.array(["a", "b", "c", "d", "e", "f", "g", "h"])
dummy_1unique = np.array([1, 1, 1, 1, 1, 1, 1, 1])
dummy_empty = np.array(
    [
        float("inf"),
        float("inf"),
        float("inf"),
        float("inf"),
        float("inf"),
        float("inf"),
        float("inf"),
        float("inf"),
    ]
)


class TestJonckheereTrendTest:
    def test_valid_result(self):
        # This test should be valid
        x = dummy_data
        group = dummy_group
        result = jonckheere(x, group, alternative="increasing")
        assert isinstance(result, dict)
        assert "statistic" in result
        assert "p.value" in result
        assert "alternative" in result
        assert result["alternative"] == "increasing"
        assert 0 <= result["p.value"] <= 1

    def test_data_not_numeric(self):
        # Data non-numeric to trigger warning
        x = dummy_nonnumeric
        group = dummy_group
        with pytest.raises(ValueError, match="Data needs to be numeric"):
            jonckheere(x, group, alternative="two.sided")

    def test_group_not_numeric(self):
        # Group non-numeric to trigger warning
        x = dummy_data
        group = dummy_nonnumeric
        with pytest.raises(ValueError, match="Group needs to be numeric or ordered factor"):
            jonckheere(x, group, alternative="two.sided")

    def test_group_data_different_lengths(self):
        # Group and data are different lengths to trigger warning
        x = dummy_data
        group = dummy_uneven
        with pytest.raises(ValueError, match="Data and group values need to be the same length"):
            jonckheere(x, group, alternative="two.sided")

    def test_alternative_not_valid(self):
        # Using an alternative value that is not valid to trigger warning
        x = dummy_data
        group = dummy_group
        alternative = "test"
        with pytest.raises(ValueError, match="Alternative choice not valid"):
            jonckheere(x, group, alternative)

    def test_data_empty(self):
        # Data is missing to trigger warning
        x = dummy_empty
        group = dummy_group
        with pytest.raises(
            ValueError, match="Either data or group is missing for all observations"
        ):
            jonckheere(x, group, alternative="two.sided")

    def test_group_empty(self):
        # Group data is missing to trigger warning
        x = dummy_data
        group = dummy_empty
        with pytest.raises(
            ValueError, match="Either data or group is missing for all observations"
        ):
            jonckheere(x, group, alternative="two.sided")

    def test_one_group(self):
        # Testing only one unique group to trigger warning
        x = dummy_data
        group = dummy_1unique
        with pytest.raises(ValueError, match="Only one group has non-missing data"):
            jonckheere(x, group, alternative="two.sided")
