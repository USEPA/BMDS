#!/bin/bash

set -xe

# only do this first time
DEPS="$GITHUB_WORKSPACE/deps"
if [ -d "$DEPS" ]; then
    exit 0
else
    mkdir -p $DEPS
fi

# eigen
cd $GITHUB_WORKSPACE
cp ./vendor/eigen-3.4.0.tar.gz $DEPS
cd $DEPS && tar -xf eigen-3.4.0.tar.gz

# nlopt
cd $GITHUB_WORKSPACE
cp ./vendor/nlopt-2.7.1.tar.gz $DEPS
cd $DEPS && tar -xf nlopt-2.7.1.tar.gz
cd nlopt-2.7.1 && mkdir build && cd build && cmake -DBUILD_SHARED_LIBS=OFF .. && cmake --build . --config Release

# gsl
cd $GITHUB_WORKSPACE
cp ./vendor/ampl-gsl-60539d2.tar.gz $DEPS
cd $DEPS && tar -xf ampl-gsl-60539d2.tar.gz
cd ampl-gsl-60539d2 && mkdir build && cd build && cmake .. -DNO_AMPL_BINDINGS=1 && cmake --build . --config Release
