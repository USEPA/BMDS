from importlib.metadata import version

from . import bleep
from .doubles import double_add

__version__ = version("pybmds")
__all__ = ["__version__", "bleep", "double_add"]
