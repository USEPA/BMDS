import numpy as np
import pytest

from pybmds.stats.normal import _ERR_MSG, exact_rnorm


class Test_exact_rnorm:
    def test_success(self):
        values, _ = exact_rnorm(10, 5, 1)
        assert np.isclose(np.mean(values), 5)
        assert np.isclose(np.std(values), 1)

    def test_failure(self):
        # cannot get a mean=0 and SD=5 w/ only positive numbers
        with pytest.raises(ValueError, match=_ERR_MSG):
            exact_rnorm(10, 0, 5, impose_positivity=True)

    def test_n_1(self):
        result = exact_rnorm(1, 10, 0)
        assert result == (10, 0)
        with pytest.raises(ValueError, match=_ERR_MSG):
            exact_rnorm(1, 10, 1)

    def test_n_2(self):
        values, iterations = exact_rnorm(2, 10, 1)
        assert np.allclose(values, [9, 11])
        assert iterations == 0
        with pytest.raises(ValueError, match=_ERR_MSG):
            exact_rnorm(2, 0, 5, impose_positivity=True)

    @pytest.mark.parametrize("tolerance", [0.001, 0.01, 0.05])
    def test_tolerance(self, tolerance):
        target_mean = 1
        target_sd = 0.1
        values, _ = exact_rnorm(10, target_mean, target_sd, tolerance=tolerance)

        calculated_mean = np.mean(values)
        calculated_sd = np.std(values)

        mean_diff = abs(calculated_mean - target_mean)
        sd_diff = abs(calculated_sd - target_sd)

        assert mean_diff < tolerance
        assert sd_diff < tolerance

    def test_positivity(self):
        # not all values will be positive with flag off
        values, _ = exact_rnorm(10, 0.5, 0.5, impose_positivity=False)
        assert not all(r > 0 for r in values)

        # but all values wll be positive with flag on
        values, _ = exact_rnorm(10, 0.5, 0.5, impose_positivity=True)
        assert all(r > 0 for r in values)
