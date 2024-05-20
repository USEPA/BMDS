# BMDS

## Disclaimer

This software/application has been approved for release by the U.S. Environmental Protection Agency (USEPA). Although the software has been subjected to rigorous review, the USEPA reserves the right to update the software as needed pursuant to further analysis and review. No warranty, expressed or implied, is made by the USEPA or the U.S. Government as to the functionality of the software and related material nor shall the fact of release constitute any such warranty. Furthermore, the software is released on condition that neither the USEPA nor the U.S. Government shall be held liable for any damages resulting from its authorized or unauthorized use.

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
