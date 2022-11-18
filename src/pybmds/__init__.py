from importlib.metadata import version

from . import bleep

__version__ = version("pybmds")
__all__ = ["__version__", "bleep", "double_add"]
