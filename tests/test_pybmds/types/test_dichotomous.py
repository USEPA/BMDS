import pytest
from pydantic import ValidationError

from pybmds.models import dichotomous
from pybmds.types.dichotomous import DichotomousModelSettings, DichotomousRiskType


class TestDichotomousAnalysisCPPStructs:
    def test_cpp_str(self, ddataset2):
        # ensure we can generate a string representation of the cpp structs
        model = dichotomous.Logistic(ddataset2)
        model.execute()
        text = str(model.structs)
        assert """- python_dichotomous_analysis""" in text
        assert """- python_dichotomous_model_result""" in text
        assert len(text.splitlines()) == 65


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
