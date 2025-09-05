import numpy as np
import pytest

from pybmds.stats.cochran_armitage import cochran_armitage


@pytest.fixture
def dataset() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    # tuple of (dose, n, incidence)
    return (
        np.array([0, 25, 75, 125, 200]),
        np.array([20, 20, 20, 20, 20]),
        np.array([0, 1, 2, 3, 5]),
    )


class TestCochranArmitage:
    def test_valid_result(self, dataset):
        dose, n, incidence = dataset
        result = cochran_armitage(dose, n, incidence)
        assert 0 <= result.p_value_asymptotic <= 1
        assert 0 <= result.p_value_exact <= 1

    def test_result_table(self, dataset):
        dose, n, incidence = dataset
        result = cochran_armitage(dose, n, incidence)
        assert isinstance(result.tbl(), str)

    def test_accurate_result(self, dataset):
        # This test determines if the result is accurate based on the given dataset
        dose, n, incidence = dataset
        result = cochran_armitage(dose, n, incidence)
        assert pytest.approx(result.statistic) == -2.7389250960176055
        assert pytest.approx(result.p_value_asymptotic) == 0.003082020765063486
        assert pytest.approx(result.p_value_exact) == 0.004414453999236329

    def test_dose_diff_length(self, dataset):
        # Dose data is a different length than N and incidences to trigger warning
        dose, n, incidence = dataset
        dose = dose[:-1]
        with pytest.raises(ValueError, match="All input vectors must be of the same length."):
            cochran_armitage(dose, n, incidence)

    def test_ns_diff_length(self, dataset):
        # N is a different length than data and incidences to trigger warning
        dose, n, incidence = dataset
        n = n[:-1]
        with pytest.raises(ValueError, match="All input vectors must be of the same length."):
            cochran_armitage(dose, n, incidence)

    def test_incidences_diff_length(self, dataset):
        # Incidences is a different length than N and data to trigger warning
        dose, n, incidence = dataset
        incidence = incidence[:-1]
        with pytest.raises(ValueError, match="All input vectors must be of the same length."):
            cochran_armitage(dose, n, incidence)

    def test_less_three_doses(self, dataset):
        # Less than three doses are used to trigger warning
        dose, n, incidence = dataset
        dose = dose[:2]
        with pytest.raises(ValueError, match="At least three dose groups are required."):
            cochran_armitage(dose, n, incidence)

    def test_dose_groups_too_small(self, dataset):
        # Dose groups have less than two individuals to trigger warning
        dose, n, incidence = dataset
        n[2] = 1
        with pytest.raises(ValueError, match="All dose groups must have at least two individuals."):
            cochran_armitage(dose, n, incidence)

    def test_incidences_negative(self, dataset):
        # Inputting a negative number of responders to trigger a warning
        dose, n, incidence = dataset
        incidence[2] = -1
        with pytest.raises(ValueError, match="The number of responders must be nonnegative."):
            cochran_armitage(dose, n, incidence)

    def test_more_incidences_than_group(self, dataset):
        # Inputting more responders than group total to trigger warning
        dose, n, incidence = dataset
        incidence[2] = n[2] + 1
        with pytest.raises(ValueError, match="The number of responders cannot exceed group total."):
            cochran_armitage(dose, n, incidence)

    def test_doses_not_increasing(self, dataset):
        # Inputting dose data that is not strictly increasing to trigger warning
        dose, n, incidence = dataset
        dose = dose[::-1]
        with pytest.raises(ValueError, match="Dose levels must be strictly increasing."):
            cochran_armitage(dose, n, incidence)
