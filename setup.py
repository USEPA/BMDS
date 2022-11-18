from glob import glob

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import find_packages, setup

__version__ = "0.0.1"

setup(
    name="pybmds",
    version=__version__,
    description="Python interface for US EPA Benchmark Dose Modeling Software",
    author="The BMDS development team",
    license="MIT",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.10",
    ],
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    cmdclass={"build_ext": build_ext},
    ext_modules=[
        Pybind11Extension(
            "pybmds.bleep",
            sorted(glob("src/pybmdscpp/*.cpp")),  # Sort for reproducibility
            define_macros=[("VERSION_INFO", __version__)],
        )
    ],
)
