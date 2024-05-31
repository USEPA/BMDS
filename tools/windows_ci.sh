#!/bin/bash

set -xe

# filenames and versions
EIGEN="eigen-3.4.0"
NLOPT="nlopt-2.7.1"
GSL="ampl-gsl-60539d2"

# only do this first time
DEPS="$GITHUB_WORKSPACE/deps"
if [ -d "$DEPS" ]; then
    exit 0
else
    mkdir -p $DEPS
fi

# eigen
cd $GITHUB_WORKSPACE
cp ./vendor/$EIGEN.tar.gz $DEPS
cd $DEPS && tar -xf $EIGEN.tar.gz && mv ./$EIGEN ./eigen

# nlopt
cd $GITHUB_WORKSPACE
cp ./vendor/$NLOPT.tar.gz $DEPS
cd $DEPS && tar -xf $NLOPT.tar.gz && mv ./$NLOPT ./nlopt
cd ./nlopt && mkdir build && cd build && cmake -DBUILD_SHARED_LIBS=OFF .. && cmake --build . --config Release

# gsl
cd $GITHUB_WORKSPACE
cp ./vendor/$GSL.tar.gz $DEPS
cd $DEPS && tar -xf $GSL.tar.gz  && mv ./$GSL ./gsl
cd ./gsl && mkdir build && cd build && cmake .. -DNO_AMPL_BINDINGS=1 && cmake --build . --config Release
