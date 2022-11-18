import pybmds


def test_bleep():
    # test calling cpp extension directly
    assert pybmds.bleep.add(2, 2) == 4

    # version -> python == cpp == expected value
    assert pybmds.bleep.__version__ == pybmds.__version__ == "0.0.1"
