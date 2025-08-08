import pytest
from pydantic import ValidationError

import pybmds
from pybmds.types.multi_tumor import MultitumorAnalysis, MultitumorSettings


class TestMultitumorSettings:
    def test_settings(self):
        MultitumorSettings(degrees=[0, 0])

    def test_no_extra(self):
        with pytest.raises(ValidationError):
            MultitumorSettings(degrees=[0, 0], foo=123)


class TestMultitumorAnalysis:
    def test_str(self, mt_datasets):
        session = pybmds.Multitumor(mt_datasets, degrees=[1, 1, 1])
        session.execute()
        assert isinstance(session.structs, MultitumorAnalysis)
        assert "python_multitumor_analysis" in str(session.structs)
