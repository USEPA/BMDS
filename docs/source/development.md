# Development

Like layers of an onion, there are multiple BMDS products that can be developed, each of which require slightly different build steps and tooling. A description of each of the environments is described below.

## Building and testing `pybmds`

To install a development version:

```bash
# clone repo
git clone https://github.com/USEPA/BMDS

# create a new python virtual environment
python -m venv venv

# active virtual environment
source venv/bin/activate

# install package in developer mode and developer tools
python -m pip install -U pip
python -m pip install -e ".[dev,docs]"
```

This will install the package and all dependencies, including building the C++ `bmdscore` shared object which will require a compiler and related dependencies installed and configured for your local setup.

:::{tip}

If you want to skip building `bmdscore` locally, you can set an environment variable `SKIP_BUILD=1` in your python environment; this will allow you to install the package but skip compilation of the shared object.

:::

Tests are written using [pytest](http://doc.pytest.org/en/latest/). To run all tests:

```bash
# run all tests
py.test

# To run a specific test
py.test -k test_my_special_test_name
```

There is a built-in Makefile, ``make``, that runs many common utility methods for developing the application. You can run tests by running `make test`, run code formatting using `make format`, or build documentation using `make docs`. For more details or other built-in actions, view the `Makefile` (or the `make.bat` file for Windows).

Github Actions are in place to execute whenever code a pull-request is submitted to check code formatting and successful tests. When code is merged into the `main` branch, a wheel artifact is created and stored on github. We use [cibuildwheel](https://cibuildwheel.pypa.io/en/stable/) to build wheel packages in Windows, Mac, and Linux using Github Actions.

## Build and testing in `bmdscore`

You can also build the `bmdscore` C++ library separately from `pybmds` integration. Inside `src/build` directory:

```bash
# build C++ code
cmake ..
make

# run all tests
make run_tests

# generate coverage report
make coverage

#run all tests and generate coverage report
make run_tests_with_coverage
```

## Build and testing `bmds-ui` (BMDS Desktop)

See [documentation](https://github.com/USEPA/BMDS-UI/blob/main/docs/development.md) in the BMDS-UI repository.

## Documentation

Documentation is generated using [Sphinx](https://www.sphinx-doc.org/). Narrative text is written in markdown; any examples of using `pybmds` are written in jupyter notebooks and are executing when building documentation; the returned notebook with any output (such as figures or output tables) are rendered directly in the documentation.

To build documentation, there are a few different commands available in the `make` file:

```bash
make docs       #  Build documentation {html}
make docs-docx  #  Build documentation {docx}
make docs-serve #  Realtime documentation preview
make docs-clean #  Clean documentation
```

Using the `make serve` command is recommended for editing documentation; it updates the preview in realtime files are saved.

### Priors Report

The `pybmds` package includes Bayesian priors and frequentist parameter initialization settings that have been tuned to help improve model fit performance. To generate a report of the settings in all permutations, run the following command:

```bash
bmds-priors-report priors_report.md
```
