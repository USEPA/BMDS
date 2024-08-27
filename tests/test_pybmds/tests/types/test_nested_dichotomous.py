import pytest
from pydantic import ValidationError

from pybmds.models.nested_dichotomous import NestedLogistic
from pybmds.types.nested_dichotomous import (
    NestedDichotomousAnalysis,
    NestedDichotomousModelSettings,
    NestedDichotomousResult,
)


class TestNestedDichotomousModelSettings:
    def test_bootstrap_seed(self):
        # check that seed is random if not specified
        settings = [NestedDichotomousModelSettings() for _ in range(5)]
        assert len({setting.bootstrap_seed for setting in settings}) > 1

        # check you can set the seed if you need to
        settings = NestedDichotomousModelSettings(bootstrap_seed=123)
        assert settings.bootstrap_seed == 123

    def test_no_extra(self):
        with pytest.raises(ValidationError):
            NestedDichotomousModelSettings(foo=123)


class TestNestedDichotomousAnalysis:
    def test_reporting(self, nd_dataset4):
        # test various reporting characteristics of the analysis
        analysis = NestedLogistic(nd_dataset4, settings=dict(bootstrap_seed=1))
        analysis.execute()
        structs = analysis.structs
        assert isinstance(structs, NestedDichotomousAnalysis)

        # check __str__ method works w/ introspection
        assert len(str(structs)) > 0

        results = analysis.results
        assert isinstance(results, NestedDichotomousResult)
        rows = results.parameter_rows(extras={"zzz": "test"})
        assert len(rows) == len(results.parameter_names)
        assert rows[0]["zzz"] == "test"
        assert rows[0]["name"] == "a"
