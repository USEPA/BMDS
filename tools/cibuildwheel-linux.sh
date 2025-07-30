#!/bin/bash

set -xe

yum update -y
yum install -y cmake gsl-devel eigen3-devel

NLOPT="nlopt-2.10.0"

cp ./vendor/$NLOPT.tar.gz ~
cd ~
tar -xf $NLOPT.tar.gz && cd $NLOPT && mkdir build && cd build && cmake -DBUILD_SHARED_LIBS=OFF .. && make install
cd $GITHUB_WORKSPACE
