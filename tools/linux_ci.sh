#!/bin/bash

sudo apt-get update -y
sudo apt-get install -y automake build-essential libtool make cmake libgslcblas0 libgsl-dev libeigen3-dev libnlopt-dev libnlopt-cxx-dev

export "EIGEN_DIR=/usr/include/eigen3"
export "NLOPT_DIR=/usr/lib/x86_64-linux-gnu/"
export "CMAKE_C_COMPILER=/usr/bin/gcc-12"
export "CMAKE_CXX_COMPILER=/usr/bin/g++-12"
export "PYTHON_EXECUTABLE=$Python3_ROOT_DIR/bin/python"
export "PYTHON_LIBRARY_DIR=$Python3_ROOT_DIR/lib/python3.11/site-packages"
