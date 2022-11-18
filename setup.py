from glob import glob

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

__version__ = "0.0.1"

setup(
    version=__version__,
    cmdclass={"build_ext": build_ext},
    ext_modules=[
        Pybind11Extension(
            "pybmds.bleep",
            sorted(glob("src/pybmdscpp/*.cpp")),  # Sort for reproducibility
            define_macros=[("VERSION_INFO", __version__)],
        )
    ],
)
