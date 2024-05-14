# BMDS

TODO - describe!

Disclaimer: The United States Government project code is provided on an "as is" basis and the user assumes responsibility for its use. The United States Government has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the United States Government. The NIH or EPA seal and logo shall not be used in any manner to imply endorsement of any commercial product or activity by NIH or EPA or the United States Government.

## ToxicR

TODO

## pybmds

### Developer setup

Make sure you have python 3.11 or higher available on your path.

```bash
# clone project
git clone git@github.com:USEPA/bmds.git
cd bmds

# create virtual environment and activate
python -m venv venv --prompt pybmds
source venv/bin/activate  # or venv\Scripts\activate on windows.

# install packages
python -m pip install -U pip
python -m pip install -r requirements_dev.txt

# build
make build     # recompile source for development of pybmds package

# test local install
pybmds hello

# code formatting
make lint      # identify formatting errors
make format    # fix formatting errors when possible

# testing
make test      # run tests
make coverage  # show test coverage

# distribution
make dist     # build a portable python wheel for distribution
```

Github actions are setup to execute whenever code is pushed to check code formatting and successful tests. In addition, when code is pushed to the `main` branch, a wheel artifact is created and stored on github.

## C++ build and testing
inside src/build directory

```bash
# build c++ code
cmake ..
make

# run all tests
make run_tests

# generate coverage report
make coverage

#run all tests and generate coverage report
make run_tests_with_coverage
```
