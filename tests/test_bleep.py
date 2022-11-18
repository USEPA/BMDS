import pybmds
from pybmds import bleep


def test_bleep():
    # test calling cpp extension directly
    assert bleep.add(2, 2) == 4
    assert bleep.sub(2, 2) == 0

    # version -> python == cpp == expected value
    assert pybmds.bleep.__version__ == pybmds.__version__ == "0.0.1"
