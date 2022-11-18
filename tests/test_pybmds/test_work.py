import pytest

from pybmds.work import bottles, double_add, super_add


def test_bottles(capsys):
    bottles(num=2, beverage="coke")
    captured = capsys.readouterr()
    assert "2 bottles of coke on the wall" in captured.out


def test_super_add():
    # test expected
    assert super_add(2, 2) == 4

    # test unexpected
    with pytest.raises(TypeError):
        super_add(2, "skunk")


def test_double_add():
    assert double_add(2, 2) == 8
