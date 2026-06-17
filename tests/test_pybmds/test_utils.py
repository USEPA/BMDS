from dataclasses import dataclass

import numpy as np

import pybmds
from pybmds.constants import BMDS_BLANK_VALUE
from pybmds.utils import ff, get_version, unique_items


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


def test_ff_formats_missing_and_nonfinite_values():
    assert ff(BMDS_BLANK_VALUE) == "-"
    assert ff(float("nan")) == "-"
    assert ff(np.inf) == "-"
