# Development

BMDS consists of multiple products, each of which requires slightly different build steps and tooling. The following sections describe each of those environments.

## Building and testing `pybmds`

For `pybmds` development, we use [uv](https://docs.astral.sh/uv/) to download Python and manage environments. Install `uv` and make sure it's available on your terminal path before proceeding.

The Python `pybmds` package's engine, `bmdscore`, is a C++ library, and therefore it must be compiled in order to be used. Compiling C++ is difficult, and may not be required for some development tasks.  Instructions below describe a few different installation paths.

If you're on Linux/MacOS, it's easier to build the C++ package, so we recommend following the guide below to build `pybmds` and compile `bmdscore` (see [below](#compiling-pybmds)).

However, on Windows, it's more complex.  The simplest method is to install the Python package without compiling C++. To do that, you'll need to download a precompiled file from GitHub and place that file in the correct location.

```ps1
# clone repository
git clone https://github.com/USEPA/BMDS

# create a new python virtual environment
uv venv --python=3.13

# activate virtual environment
.venv/Scripts/activate

# download a Python `bmdscore` artifact from GitHub:
# https://github.com/USEPA/BMDS/actions/workflows/test-windows.yml
# move the `*.pyd` file to the `src/pybmds` folder, beside the `bmdscore.pyi` file.

# install pybmds in developer mode (skipping compilation)
$env:SKIP_BUILD=1
uv pip install -e ".[dev,docs]"

# tests should pass
pytest
```
(compiling-pybmds)=
### Building `pybmds` including compiling `bmdscore`

The [vcpkg](https://github.com/microsoft/vcpkg) C++ package manager is used to install C++ dependencies. The steps below will install `vcpkg`, download and compile the required dependencies as static libraries, and then build the `bmdscore` package and Python bindings. These steps are fully automated on GitHub using GitHub Actions, so should you have any troubles for any particular environment, refer to the GitHub Actions.

For **Linux/MacOS**:

```bash
# clone repository
git clone https://github.com/USEPA/BMDS
cd bmds

# check out a snapshot from vcpkg git repository
mkdir vcpkg
cd vcpkg
git init
git remote add origin git@github.com:microsoft/vcpkg.git
git fetch --depth 1 origin c9c17dcea3016bc241df0422e82b8aea212dcb93
git checkout FETCH_HEAD
cd ..

# build static dependencies (eigen, nlopt, gsl)
# use the appropriate VCPKG_HOST_TRIPLET for your architecture/OS
export VCPKG_HOST_TRIPLET="arm64-osx" # for macOS
export VCPKG_HOST_TRIPLET="x64-linux" # for linux
./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install --host-triplet="$VCPKG_HOST_TRIPLET"

# install pybind11 at the root so for our build phase
uv pip install pybind11=="3.0.0" --target="./pybind11"

# set environment variables for building python extension
export CMAKE_PREFIX_PATH="$(readlink -f ./pybind11/pybind11/share/cmake)" --overlay-ports="./vendor/ports"
export CMAKE_BUILD_PARALLEL_LEVEL=$(nproc)

echo "$VCPKG_HOST_TRIPLET"
echo "$CMAKE_PREFIX_PATH"
echo "$CMAKE_BUILD_PARALLEL_LEVEL"

# create a new python virtual environment
uv venv --python="3.13"
source .venv/bin/activate

# compile and install the package
uv pip install -v -e ".[dev]"

# test to see if compilation was successful
pytest
```

For **Windows (PowerShell)**:

Install the command line tools for [Visual Studio](https://visualstudio.microsoft.com/downloads/) (not [Visual Studio Code](https://code.visualstudio.com/)). In the installation process, install the C++ dependencies that are needed to build C++ code (cmake, make, etc.). Using [Windows Terminal](https://learn.microsoft.com/en-us/windows/terminal/), select a Visual Studio Window so the environment is configured for compiling C++.

```ps1
# clone repository
git clone https://github.com/USEPA/BMDS
cd bmds

# check out a single commit from a git repository
mkdir vcpkg
cd vcpkg
git init
git remote add origin git@github.com:microsoft/vcpkg.git
git fetch --depth 1 origin c9c17dcea3016bc241df0422e82b8aea212dcb93
git checkout FETCH_HEAD
cd ..

# build static dependencies (eigen, nlopt, gsl)
$env:VCPKG_HOST_TRIPLET="x64-windows-static-md"
.\vcpkg\bootstrap-vcpkg.bat
.\vcpkg\vcpkg install --host-triplet="$env:VCPKG_HOST_TRIPLET" --overlay-ports="./vendor/ports"

# install pybind11 at the root so for our build phase
uv pip install pybind11=="3.0.0" --target="./pybind11"

# set environment variables for building python extension
$env:CMAKE_PREFIX_PATH=Resolve-Path "./pybind11/pybind11/share/cmake" | Select-Object -ExpandProperty Path
$env:CMAKE_BUILD_PARALLEL_LEVEL=[Environment]::ProcessorCount

echo "$env:VCPKG_HOST_TRIPLET"
echo "$env:CMAKE_PREFIX_PATH"
echo "$env:CMAKE_BUILD_PARALLEL_LEVEL"

# create a new python virtual environment
uv venv --python=3.13
.venv/Scripts/activate

# compile and install the package
uv pip install -v -e ".[dev]"

# test to see if compilation was successful
pytest
```

## Testing `pybmds`

Tests are written using [pytest](http://doc.pytest.org/en/latest/). To run all tests:

```bash
# run all tests
py.test

# To run a specific test
py.test -k test_my_special_test_name
```

You can run tests by running `poe test`, run code formatting using `poe format`, or build documentation using `poe docs`. For more details on other built-in actions, view the `poe` commands by simply typing `poe`.

Github Actions are in place to execute whenever a pull-request is submitted to check code formatting and successful tests. When code is merged into the `main` branch, a wheel artifact is created and stored on github. We use [cibuildwheel](https://cibuildwheel.pypa.io/en/stable/) to build wheel packages in Windows, macOS, and Linux using Github Actions.

## Build and testing `bmdscore`

You can also build the `bmdscore` C++ library separately from `pybmds` integration. Inside `src/build` directory:

```bash
# check out a single commit from a git repository
mkdir vcpkg
cd vcpkg
git init
git remote add origin git@github.com:microsoft/vcpkg.git
git fetch --depth 1 origin c9c17dcea3016bc241df0422e82b8aea212dcb93
git checkout FETCH_HEAD
cd ..

# build static dependencies (eigen, nlopt, gsl)
$env:VCPKG_HOST_TRIPLET="x64-linux". # OS-specific; x64-windows-static-md / arm64-osx / etc.
.\vcpkg\bootstrap-vcpkg
.\vcpkg\vcpkg install --host-triplet="$env:VCPKG_HOST_TRIPLET" --overlay-ports="./vendor/ports"

# build C++ code
mkdir -p src/build
cd src/build
cmake ..
make

# run all tests
make run_tests

# generate coverage report
make coverage

#run all tests and generate coverage report
make run_tests_with_coverage
```

Code formatting is also enforced in the C++ codebase. To run the formatter, use the following command:

```bash
poe format-cpp
```

A specific version of clang-format is pinned for reproducibility on all platforms.

## Build and testing `bmds-ui` (BMDS Desktop)

See [documentation](https://github.com/USEPA/BMDS-UI/blob/main/docs/development.md) in the BMDS-UI repository.

## Documentation

Documentation is generated using [Sphinx](https://www.sphinx-doc.org/). Narrative text is written in Markdown; any examples of using `pybmds` are written in jupyter notebooks and are executing when building documentation; the returned notebook with any output (such as figures or output tables) are rendered directly in the documentation.

To build documentation, there are a few different commands available:

```bash
poe docs       #  Build documentation {html}
poe docs-docx  #  Build documentation {docx}
poe docs-serve #  Realtime documentation preview
```

Using the `poe docs-serve` command is recommended for editing documentation; it updates the preview in realtime files are saved.

## Versioning

We use [calendar versioning](https://calver.org/) for `pybmds`, where:

* `major` is the year of the release (ex: `28` for a 2028 release)
* `minor` is incremented for each release of the calendar year, starting at `1`
* `aN` is the alpha release for testing, where N starts at `1`
* `dev` is any upcoming pre-release currently under development.

As an example, consider the scenario where we're beginning development our first release in 2028:

* In `pybmds.__init__`, set `__version__ = "28.1a1.dev"`
* Iterate until we're ready for an alpha release
    * Update the version to `28.1a1`, and git tag the release `28.1a1`
    * Immediately change the `main` branch to `28.1a2.dev`
* Begin testing of `28.1a1`
    * If changes are needed, iterate on `28.1a2.dev`
    * If changes are not needed, release a `28.1` by changing the version and minting a tag

The [packaging](https://packaging.pypa.io/en/stable/index.html) package implements [PEP440](https://peps.python.org/pep-0440/), and can be used to check candidate versions:

```python
from packaging.version import Version

Version('28.1a1.dev')
# _Version(release=(28, 1), pre=('a', 1), dev=('dev', 0))
```

### Priors Report

The `pybmds` package includes Bayesian priors and frequentist parameter initialization settings that have been tuned to help improve model fit performance. To generate a report of the settings in all permutations, run the following command:

```bash
bmds-priors-report priors_report.md
```
