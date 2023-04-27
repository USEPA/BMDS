#!/bin/bash

sudo apt-get update -y
sudo apt-get install -y automake build-essential libtool make cmake libgslcblas0 libgsl-dev libeigen3-dev libnlopt-dev

echo "a"
/usr/bin/gsl-config --libs
echo "b"
/usr/bin/gsl-config --prefix
echo "c"
/usr/bin/gsl-config --libs-without-cblas

export "EIGEN_DIR=/usr/include/eigen3/Eigen/"
export "GSL_DIR=/usr/include/gsl"
export "NLOPT_DIR=/usr/include"
export "CMAKE_C_COMPILER=/usr/bin/gcc-12"
export "CMAKE_CXX_COMPILER=/usr/bin/g++-12"
export "PYTHON_EXECUTABLE=/opt/hostedtoolcache/Python/3.11.3/x64/bin/python"
export "PYTHON_LIBRARY_DIR=/opt/hostedtoolcache/Python/3.11.3/x64/lib/python3.11/site-packages"
