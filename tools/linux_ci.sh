#!/bin/bash

sudo apt-get update -y
sudo apt-get install -y automake build-essential libtool make cmake libgslcblas0 libgsl-dev libeigen3-dev libnlopt-dev

export "EIGEN_DIR=/usr/include/eigen3"
export "NLOPT_DIR=/usr/lib/x86_64-linux-gnu/"
export "CMAKE_C_COMPILER=/usr/bin/gcc-12"
export "CMAKE_CXX_COMPILER=/usr/bin/g++-12"
export "PYTHON_EXECUTABLE=/opt/hostedtoolcache/Python/3.11.3/x64/bin/python"
export "PYTHON_LIBRARY_DIR=/opt/hostedtoolcache/Python/3.11.3/x64/lib/python3.11/site-packages"
