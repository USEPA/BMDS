# Development

BMDS consists of multiple products, each of which requires slightly different build steps and tooling. The following sections describe each of those environments.

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

This will install the package and all dependencies, including building the C++ `bmdscore` shared object, which will require a compiler and related dependencies to be installed and configured for your local setup.

:::{tip}

If you want to skip building `bmdscore` locally, you can set the environment variable `SKIP_BUILD=1` in your Python environment. Setting this variable will allow you to install the package but skip compilation of the shared object.

:::

Tests are written using [pytest](http://doc.pytest.org/en/latest/). To run all tests:

```bash
# run all tests
py.test

# To run a specific test
py.test -k test_my_special_test_name
```

You can run tests by running `poe test`, run code formatting using `poe format`, or build documentation using `poe docs`. For more details on other built-in actions, view the `poe` commands by simply typing `poe`.

Github Actions are in place to execute whenever a pull-request is submitted to check code formatting and successful tests. When code is merged into the `main` branch, a wheel artifact is created and stored on github. We use [cibuildwheel](https://cibuildwheel.pypa.io/en/stable/) to build wheel packages in Windows, macOS, and Linux using Github Actions.

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
