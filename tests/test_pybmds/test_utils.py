from dataclasses import dataclass

import pybmds
from pybmds.utils import get_version, unique_items


def test_citation():
    assert "pybmds" in pybmds.citation()


def test_get_version():
    version = get_version()
    assert int(version.dll.split(".")[0]) >= 24  # assume dll in format "YY.MM..."


@dataclass
class Foo:
    bar: str


def test_unique_items():
    assert unique_items([Foo(bar="b"), Foo(bar="a"), Foo(bar="b")], "bar") == "a, b"
