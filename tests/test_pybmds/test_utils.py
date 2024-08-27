import pybmds
from pybmds.utils import get_version


def test_citation():
    assert pybmds.citation().startswith("pybmds.")


def test_get_version():
    version = get_version()
    assert int(version.dll.split(".")[0]) >= 24  # assume dll in format "YY.MM..."
