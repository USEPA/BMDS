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

Tests are written using [pytest](http://doc.pytest.org/en/latest/). To run all tests:

```bash
# run all tests
py.test

# To run a specific test
py.test -k test_my_special_test_name
```

There is a built in Makefile command, ``make dev`` that creates a tmux application which auto-update the documentation; check out the ``Makefile`` for a list of other built-in actions.

## Priors Report

The pybmds package includes Bayesian priors and frequentist parameter initialization settings that have been tuned to help improve model fit performance. To generate a report of the settings in all the possible permutations, run the command:

```bash
bmds-priors-report priors_report.md
```
