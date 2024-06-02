#!/bin/bash

set -xe

# only do this first time
DEPS="$GITHUB_WORKSPACE/deps"
if [ -d "$DEPS" ]; then
    exit 0
fi
mkdir -p $DEPS

# filenames and versions
EIGEN="eigen-3.4.0"
NLOPT="nlopt-2.7.1"
GSL="ampl-gsl-60539d2"

# eigen
cp $GITHUB_WORKSPACE/vendor/$EIGEN.tar.gz $DEPS
cd $DEPS && tar -xf $EIGEN.tar.gz && mv $EIGEN eigen && rm $EIGEN.tar.gz

# nlopt
cp $GITHUB_WORKSPACE/vendor/$NLOPT.tar.gz $DEPS
cd $DEPS && tar -xf $NLOPT.tar.gz && mv $NLOPT nlopt && rm $NLOPT.tar.gz
cd nlopt && mkdir build && cd build
cmake -DBUILD_SHARED_LIBS=OFF ..
cmake --build . --config Release

# gsl
cp $GITHUB_WORKSPACE/vendor/$GSL.tar.gz $DEPS
cd $DEPS && tar -xf $GSL.tar.gz && mv gsl-master gsl && rm $GSL.tar.gz
cd gsl && mkdir build && cd build
cmake -DNO_AMPL_BINDINGS=1 ..
cmake --build . --config Release
