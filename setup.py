from glob import glob

from pybind11.setup_helpers import ParallelCompile, Pybind11Extension, build_ext, naive_recompile
from setuptools import find_packages, setup

__version__ = "0.0.1"

ParallelCompile("NPY_NUM_BUILD_JOBS", needs_recompile=naive_recompile).install()

setup(
    name="pybmds",
    version=__version__,
    cmdclass={"build_ext": build_ext},
    license="MIT",
    install_requires=[
        "typer",
    ],
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "pybmds = pybmds.cli:app",
        ]
    },
    ext_modules=[
        Pybind11Extension(
            "pybmds.bleep",
            sorted(glob("src/pybmdscpp/*.cpp")),  # Sort for reproducibility
            define_macros=[("VERSION_INFO", __version__)],
        )
    ],
)
