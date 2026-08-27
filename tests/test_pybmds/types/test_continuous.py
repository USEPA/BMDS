import warnings

import numpy as np
import pytest
from pydantic import ValidationError

from pybmds.constants import DistType, PriorClass
from pybmds.models import continuous
from pybmds.types.continuous import (
    ContinuousModelSettings,
    ContinuousParameters,
    ContinuousRiskType,
    _display_blank_value,
    _display_loud_loglikelihood,
    _display_loud_summary_value,
)


def test_loud_display_helpers_cover_numeric_and_text_values():
    assert _display_blank_value(np.nan) == "-"
    assert _display_blank_value("value") == "value"
    assert _display_loud_summary_value(1.25) == "1.25"
    assert _display_loud_summary_value("value") == "value"
    assert _display_loud_loglikelihood(1.25) == "-1.25"
    assert _display_loud_loglikelihood("value") == "value"


def test_continuous_parameters_loud_draw_shape_validation():
    model = type(
        "Model",
        (),
        {
            "get_param_names": lambda self: ["a", "b"],
            "get_priors_list": lambda self: [[0, 1, 1, -1, 1], [0, 1, 1, -1, 1]],
        },
    )()
    assert ContinuousParameters.from_loud_draws(model, np.array([1.0, 2.0])).values.shape == (1,)
    assert ContinuousParameters.from_loud_draws(model, np.array([[1.0, 2.0]])).cov.shape == (2, 2)
    with pytest.raises(ValueError, match="Unsupported LOUD parameter draw shape"):
        ContinuousParameters.from_loud_draws(model, np.zeros((1, 1, 1)))
    with pytest.raises(ValueError, match="draws are empty"):
        ContinuousParameters.from_loud_draws(model, np.empty((1, 0)))


class TestContinuousModelSettings:
    @pytest.mark.parametrize(
        "bmr,bmr_type,expected_text",
        [
            (1.0, ContinuousRiskType.StandardDeviation, "1.0 Standard Deviation"),
            (0.1, ContinuousRiskType.RelativeDeviation, "10% Relative Deviation"),
            (1, ContinuousRiskType.AbsoluteDeviation, "1.0 Absolute Deviation"),
            (1, ContinuousRiskType.HybridAdded, "1.0 Hybrid Added"),
            (1, ContinuousRiskType.HybridExtra, "1.0 Hybrid Extra"),
            (1, ContinuousRiskType.Extra, "1.0 Extra"),
            (1, ContinuousRiskType.HybridExtra, "1.0 Hybrid Extra"),
            (1, ContinuousRiskType.HybridAdded, "1.0 Hybrid Added"),
        ],
    )
    def test_bmr_text(self, bmr, bmr_type, expected_text):
        settings = ContinuousModelSettings(bmr=bmr, bmr_type=bmr_type)
        assert settings.bmr_text == expected_text

    def test_no_extra(self):
        with pytest.raises(ValidationError):
            ContinuousModelSettings(foo=123)

    def test_loud_mcmc_settings(self):
        settings = ContinuousModelSettings(n_chains=4, seed=123)
        assert settings.n_chains == 4
        assert settings.seed == 123

        with pytest.raises(ValidationError):
            ContinuousModelSettings(n_chains=0)

        with pytest.raises(ValidationError):
            ContinuousModelSettings(seed=-1)

    def test_loud_mcmc_sample_limits(self):
        settings = ContinuousModelSettings(
            priors=PriorClass.bayesian_loud,
            samples=25_000,
            burnin=5_000,
            n_chains=4,
        )
        assert settings.samples == 25_000
        assert settings.burnin == 5_000
        assert settings.n_chains == 4

        with pytest.raises(ValidationError, match="less than or equal to 100000"):
            ContinuousModelSettings(priors=PriorClass.bayesian_loud, samples=100_001)

        with pytest.raises(ValidationError, match="total samples across all chains"):
            ContinuousModelSettings(
                priors=PriorClass.bayesian_loud,
                samples=25_001,
                n_chains=4,
            )

        with pytest.raises(ValidationError, match="LOUD burnin cannot exceed 20%"):
            ContinuousModelSettings(
                priors=PriorClass.bayesian_loud,
                samples=10_000,
                burnin=2_001,
            )


class TestContinuousGof:
    def test_collapse(self, cdataset, cidataset):
        # goodness of fit should collapse into non-zero fields

        # continuous summary data already collapsed; no change in length of table
        model = continuous.Power(cdataset)
        res = model.execute()
        assert res.gof.n() == cdataset.num_dose_groups == 5
        assert res.gof.n() == len(cdataset.doses) == 5

        # continuous individual summary data collapses appropriately
        model = continuous.Power(cidataset)
        res = model.execute()
        assert res.gof.n() == cidataset.num_dose_groups == 7
        assert res.gof.n() == len(cidataset.doses)
        assert res.gof.n() < len(cidataset.individual_doses)
        assert res.gof.n() == len(set(cidataset.individual_doses))


class TestContinuousParameters:
    def test_from_loud_draws_uses_joint_median_draw(self, cdataset):
        model = continuous.Power(
            cdataset,
            settings=dict(disttype=DistType.normal, priors=PriorClass.bayesian_loud),
        )
        draws = np.array(
            [
                [0.0, 0.0, 0.0, 0.0],
                [10.0, 10.0, 10.0, 10.0],
                [100.0, 1.0, 100.0, 1.0],
            ]
        )

        params = ContinuousParameters.from_loud_draws(model, draws)

        assert any(np.array_equal(params.values, draw) for draw in draws)
        assert not np.allclose(params.values, np.mean(draws, axis=0))
        assert not np.allclose(params.values, np.median(draws, axis=0))

    def test_from_loud_draws_suppresses_extreme_draw_warnings(self, cdataset):
        model = continuous.Power(
            cdataset,
            settings=dict(disttype=DistType.normal, priors=PriorClass.bayesian_loud),
        )
        draws = np.array(
            [
                [1e200, 1e200, 1e200, 1e200],
                [1.1e200, 1.1e200, 1.1e200, 1.1e200],
                [1.2e200, 1.2e200, 1.2e200, 1.2e200],
            ]
        )

        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always", RuntimeWarning)
            params = ContinuousParameters.from_loud_draws(model, draws)

        assert len(params.values) == 4
        assert not [warning for warning in caught if issubclass(warning.category, RuntimeWarning)]

    def test_exp3(self, cdataset):
        """
        Edge case for exp3 - the dll expects a prior for the c parameter, but the
        returned output effectively drops the c array and shifts all other values down one.
        We check that the input and output values are shifted as required.
        """
        model = continuous.ExponentialM3(cdataset, settings=dict(disttype=DistType.normal))
        res = model.execute()
        # param names for prior are as expected
        assert model.get_param_names() == ["a", "b", "c", "d", "alpha"]
        # but outputs have been shifted
        assert res.parameters.names == ["a", "b", "d", "alpha"]

        model = continuous.ExponentialM3(cdataset, settings=dict(disttype=DistType.normal_ncv))
        res = model.execute()
        # param names for prior are as expected
        assert model.get_param_names() == ["a", "b", "c", "d", "rho", "alpha"]
        # but outputs have been shifted
        assert res.parameters.names == ["a", "b", "d", "rho", "alpha"]

        # confirm arrays all the same length after changes
        params = res.parameters
        n_params = len(params.names)
        for field in [
            "values",
            "se",
            "bounded",
            "prior_type",
            "prior_initial_value",
            "prior_stdev",
            "prior_min_value",
            "prior_max_value",
        ]:
            assert getattr(params, field).size == n_params
        assert params.cov.shape == (n_params, n_params)


class TestContinuousAnalysisCPPStructs:
    def test_cpp_str(self, cdataset2):
        # ensure we can generate a string representation of the cpp structs
        model = continuous.Power(cdataset2)
        model.execute()
        text = str(model.structs)
        assert "python_continuous_analysis" in text
        assert "python_continuous_model_result" in text
