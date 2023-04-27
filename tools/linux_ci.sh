#!/bin/bash

sudo apt-get update -y
sudo apt-get install -y automake build-essential libtool make cmake libgslcblas0 libgsl-dev libeigen3-dev libnlopt-dev

which python

ls /usr/include/eigen3/Eigen/
export "EIGEN_DIR=/usr/include/eigen3/Eigen/"
ls /usr/include/gsl
export "GSL_DIR=/usr/include/gsl"
ls /usr/include
export "NLOPT_DIR=/usr/include"
ls /usr/bin
export "CMAKE_C_COMPILER=/usr/bin/gcc"
export "CMAKE_CXX_COMPILER=/usr/bin/g++"
ls /opt/hostedtoolcache/Python/3.11.3/x64/bin
export "PYTHON_EXECUTABLE=/opt/hostedtoolcache/Python/3.11.3/x64/bin/python"

ls /opt/hostedtoolcache/Python/3.11.3/
ls /opt/hostedtoolcache/Python/3.11.3/x64
ls /opt/hostedtoolcache/Python/3.11.3/lib
ls /opt/hostedtoolcache/Python/3.11.3/lib/site-packages

export "PYTHON_LIBRARY_DIR=/opt/hostedtoolcache/Python/3.11.3/x64/lib/site-packages
