import math
import re
from collections.abc import Iterable
from dataclasses import dataclass
from datetime import datetime
from functools import cache

import numpy as np
import tabulate

from .constants import BMDS_BLANK_VALUE


def multi_lstrip(txt: str) -> str:
    """Left-strip all lines in a multiline string."""
    return "\n".join(line.lstrip() for line in txt.splitlines()).strip()


def pretty_table(data, headers):
    return tabulate.tabulate(data, headers=headers, tablefmt="fancy_outline")


def ff(value) -> str:
    """Float formatter for floats and float-like values"""
    if isinstance(value, str):
        return value
    elif abs(value) > 1e6:
        return f"{value:.1E}"
    elif value > 0 and value < 0.001:
        return "<0.001"
    elif np.isclose(value, int(value)):
        return str(int(value))
    else:
        return f"{value:.3f}".rstrip("0")


def four_decimal_formatter(value: float) -> str:
    # Expected values between 0 and 100; with happy case, returns 4 decimals
    if value == BMDS_BLANK_VALUE or not math.isfinite(value):
        return "-"
    elif value == 0:
        return "0"
    elif abs(value) > 100:
        return ff(value)
    elif value > 0 and value < 0.0001:
        return "<0.0001"
    else:
        return f"{value:.4f}".rstrip("0")


def str_list(items: Iterable) -> str:
    return ",".join([str(item) for item in items])


def citation() -> str:
    """
    Return a citation for the software.
    """
    executed = datetime.now().strftime("%B %d, %Y")
    version = get_version()
    url = "https://pypi.org/project/pybmds/"
    return f"pybmds. (Version {version.python}; Model Library Version {version.dll}) [Python package]. Available from {url}. Executed on {executed}."


@dataclass
class Version:
    python: str
    dll: str


@cache
def get_version():
    from . import __version__, bmdscore

    return Version(dll=bmdscore.version(), python=__version__)


def camel_to_title(txt: str) -> str:
    return re.sub(r"(?<=\w)([A-Z])", r" \1", txt)
