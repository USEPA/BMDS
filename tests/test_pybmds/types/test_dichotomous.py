import pytest
from pydantic import ValidationError

from pybmds.constants import PriorClass
from pybmds.models import dichotomous
from pybmds.types.dichotomous import (
    DichotomousModelSettings,
    DichotomousRiskType,
    _display_loud_loglikelihood,
    _display_loud_summary_value,
)


def test_loud_display_helpers_cover_numeric_blank_and_text_values():
    assert _display_loud_summary_value(float("nan")) == "-"
    assert _display_loud_summary_value(1.25) == "1.25"
    assert _display_loud_summary_value("value") == "value"
    assert _display_loud_loglikelihood(float("inf")) == "-"
    assert _display_loud_loglikelihood(1.25) == "-1.25"
    assert _display_loud_loglikelihood("value") == "value"


class TestDichotomousAnalysisCPPStructs:
    def test_cpp_str(self, ddataset2):
        # ensure we can generate a string representation of the cpp structs
        model = dichotomous.Logistic(ddataset2)
        model.execute()
        text = str(model.structs)
        assert "python_dichotomous_analysis" in text
        assert "python_dichotomous_model_result" in text


class TestDichotomousModelSettings:
    @pytest.mark.parametrize(
        "bmr,bmr_type,expected_text",
        [
            (0.1, DichotomousRiskType.ExtraRisk, "10% Extra Risk"),
            (0.05, DichotomousRiskType.AddedRisk, "5% Added Risk"),
        ],
    )
    def test_bmr_text(self, bmr, bmr_type, expected_text):
        settings = DichotomousModelSettings(bmr=bmr, bmr_type=bmr_type)
        assert settings.bmr_text == expected_text

    def test_no_extra(self):
        with pytest.raises(ValidationError):
            DichotomousModelSettings(foo=123)

    def test_loud_mcmc_settings(self):
        settings = DichotomousModelSettings(n_chains=4, seed=123)
        assert settings.n_chains == 4
        assert settings.seed == 123

        with pytest.raises(ValidationError):
            DichotomousModelSettings(n_chains=0)

        with pytest.raises(ValidationError):
            DichotomousModelSettings(seed=-1)

    def test_loud_mcmc_sample_limits(self):
        settings = DichotomousModelSettings(
            priors=PriorClass.bayesian_loud,
            samples=25_000,
            burnin=5_000,
            n_chains=4,
        )
        assert settings.samples == 25_000
        assert settings.burnin == 5_000
        assert settings.n_chains == 4

        with pytest.raises(ValidationError, match="less than or equal to 100000"):
            DichotomousModelSettings(priors=PriorClass.bayesian_loud, samples=100_001)

        with pytest.raises(ValidationError, match="total samples across all chains"):
            DichotomousModelSettings(
                priors=PriorClass.bayesian_loud,
                samples=25_001,
                n_chains=4,
            )

        with pytest.raises(ValidationError, match="LOUD burnin cannot exceed 20%"):
            DichotomousModelSettings(
                priors=PriorClass.bayesian_loud,
                samples=10_000,
                burnin=2_001,
            )
