# Development

To install a development version:

```bash
# clone repo
git clone https://github.com/USEPA/BMDS

# create a new python virtual environment
python -m venv venv

# active virtual environment
source venv/bin/activate

# install package in developer mode and developer tools
pip install -U pip uv
uv pip install -e ".[dev,docs]"
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

There is a built-in Makefile, ``make``, that runs many common utility methods for developing the application. For example, you can run tests by running `make test`, run code formatting using `make format`, or build documentation using `make docs`. For more details or other built-in actions, view the `Makefile` (or the `make.bat` file for Windows).

## Priors Report

The `pybmds` package includes Bayesian priors and frequentist parameter initialization settings that have been tuned to help improve model fit performance. To generate a report of the settings in all permutations, run the following command:

```bash
bmds-priors-report priors_report.md
```
