from pybmds import bmdscore, __version__


def test_version_consistency():
    # check that python and bmdscore versions are the expected value
    assert __version__ == "2023.10a1"
    assert bmdscore.version() == "2023.10a1"


def test_rbmds():
    # test calling cpp extension directly
    assert bmdscore.add2(2, 2) == 4
