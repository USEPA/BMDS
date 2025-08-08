import numpy as np
import pytest

from pybmds.stats.cochran_armitage import cochran_armitage

dummy_dose = np.array([0, 25, 75, 125, 200])
dummy_N = np.array([20, 20, 20, 20, 20])
dummy_incidence = np.array([0, 1, 2, 3, 5])
dummy_dose_short = np.array([0, 25, 75, 125])
dummy_N_short = np.array([20, 20, 20, 20])
dummy_incidence_short = np.array([0, 1, 2, 3])
dummy_dose_two = np.array([0, 25])
dummy_N_one = np.array([20, 20, 1, 20, 20])
dummy_incidence_neg = np.array([0, 1, -2, 3, 5])
dummy_incidence_large = np.array([0, 1, 2, 3, 25])
dummy_dose_decreasing = np.array([0, 25, 75, 200, 125])


class TestCochranArmitageTrendTest:
    def test_valid_result(self):
        # This test should be valid
        dose = dummy_dose
        N = dummy_N
        incidence = dummy_incidence
        result = cochran_armitage(dose, N, incidence)
        assert isinstance(result, dict)
        assert "test_statistic" in result
        assert "asymptotic_pvalue" in result
        assert "exact_pvalue" in result
        assert 0 <= result["asymptotic_pvalue"] <= 1
        assert 0 <= result["exact_pvalue"] <= 1

    def test_accurate_result(self):
        # This test determines if the result is accurate based on the given dataset
        dose = dummy_dose
        N = dummy_N
        incidence = dummy_incidence
        result = cochran_armitage(dose, N, incidence)
        assert result["test_statistic"] == -2.7389250960176055
        assert result["asymptotic_pvalue"] == 0.003082020765063486
        assert result["exact_pvalue"] == 0.004414453999236329

    def test_dose_diff_length(self):
        # Dose data is a different length than N and incidences to trigger warning
        dose = dummy_dose_short
        N = dummy_N
        incidence = dummy_incidence
        with pytest.raises(ValueError, match="All input vectors must be of the same length."):
            cochran_armitage(dose, N, incidence)

    def test_ns_diff_length(self):
        # N is a different length than data and incidences to trigger warning
        dose = dummy_dose
        N = dummy_N_short
        incidence = dummy_incidence
        with pytest.raises(ValueError, match="All input vectors must be of the same length."):
            cochran_armitage(dose, N, incidence)

    def test_incidences_diff_length(self):
        # Incidences is a different length than N and data to trigger warning
        dose = dummy_dose
        N = dummy_N
        incidence = dummy_incidence_short
        with pytest.raises(ValueError, match="All input vectors must be of the same length."):
            cochran_armitage(dose, N, incidence)

    def test_less_three_doses(self):
        # Less than three doses are used to trigger warning
        dose = dummy_dose_two
        N = dummy_N
        incidence = dummy_incidence
        with pytest.raises(ValueError, match="At least three dose groups are required."):
            cochran_armitage(dose, N, incidence)

    def test_dose_groups_too_small(self):
        # Dose groups have less than two individuals to trigger warnining
        dose = dummy_dose
        N = dummy_N_one
        incidence = dummy_incidence
        with pytest.raises(ValueError, match="All dose groups must have at least two individuals."):
            cochran_armitage(dose, N, incidence)

    def test_incidences_negative(self):
        # Inputting a negative number of responders to trigger a warning
        dose = dummy_dose
        N = dummy_N
        incidence = dummy_incidence_neg
        with pytest.raises(ValueError, match="The number of responders must be nonnegative."):
            cochran_armitage(dose, N, incidence)

    def test_more_incidences_than_group(self):
        # Inputting more responders than group total to trigger warning
        dose = dummy_dose
        N = dummy_N
        incidence = dummy_incidence_large
        with pytest.raises(ValueError, match="The number of responders cannot exceed group total."):
            cochran_armitage(dose, N, incidence)

    def test_doses_not_increasing(self):
        # Inputting dose data that is not strictly increasing to trigger warning
        dose = dummy_dose_decreasing
        N = dummy_N
        incidence = dummy_incidence
        with pytest.raises(ValueError, match="Dose levels must be strictly increasing."):
            cochran_armitage(dose, N, incidence)
