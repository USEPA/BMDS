from pybmds import RBMDS, __version__


def test_version_consistency():
    # check that python and bmdscore versions are the expected value
    assert __version__ == RBMDS.version() == "2023.10a1"


def test_rbmds():
    # test calling cpp extension directly
    assert RBMDS.add2(2, 2) == 4
