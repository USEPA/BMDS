import pytest
from pydantic import ValidationError

from pybmds.types.nested_dichotomous import NestedDichotomousModelSettings


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
