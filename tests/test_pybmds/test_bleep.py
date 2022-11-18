import pybmds


def test_pybmds():
    # test calling cpp extension directly
    assert pybmds.bleep.add(2, 2) == 4
    assert pybmds.bleep.sub(2, 2) == 0

    # test python package around c extension
    assert pybmds.double_add(2, 2) == 8

    # version -> python == cpp == expected value
    assert pybmds.bleep.__version__ == pybmds.__version__ == "0.0.1"
