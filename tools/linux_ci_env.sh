#!/bin/bash

export EIGEN_DIR="/usr/include/eigen3"
export NLOPT_DIR="/usr/lib/x86_64-linux-gnu/"
export CMAKE_C_COMPILER="/usr/bin/gcc-11"
export CMAKE_CXX_COMPILER="/usr/bin/g++-11"
export PYTHON_EXECUTABLE=$(which python)
export PYTHON_LIBRARY_DIR=$(python -c "import site; print(site.getsitepackages()[0])")
export TEST_EXC="'/usr/include/*'"
