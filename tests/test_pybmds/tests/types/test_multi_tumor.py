import pytest
from pydantic import ValidationError

from pybmds.types.multi_tumor import MultitumorSettings


class TestMultitumorSettings:
    def test_settings(self):
        MultitumorSettings(degrees=[0, 0])

    def test_no_extra(self):
        with pytest.raises(ValidationError):
            MultitumorSettings(degrees=[0, 0], foo=123)
