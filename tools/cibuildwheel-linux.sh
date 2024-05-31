#!/bin/bash

set -xe

yum update -y
yum install -y cmake gsl-devel eigen3-devel

cp ./vendor/nlopt-2.7.1.tar.gz ~
cd ~
tar -xf nlopt-2.7.1.tar.gz && cd nlopt-2.7.1 && mkdir build && cd build && cmake -DBUILD_SHARED_LIBS=OFF .. && make install
cd $GITHUB_WORKSPACE