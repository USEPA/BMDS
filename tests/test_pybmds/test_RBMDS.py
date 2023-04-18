from pybmds import RBMDS 


def test_rbmds():
    # test calling cpp extension directly
    assert RBMDS.add2(2, 2) == 4

    # version -> python == cpp == expected value
    assert RBMDS.version()  == "2023.03.1"
