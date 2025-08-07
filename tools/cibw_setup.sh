#!/bin/bash

# exit if an error occurs
set -xe

env

pwd

ls -laht /project 2>/dev/null
ls -laht /project/vcpkg_installed/ 2>/dev/null
ls -laht /project/vcpkg_installed/x64-linux-dynamic 2>/dev/null
ls -laht /project/vcpkg_installed/x64-linux-dynamic/lib 2>/dev/null
ls -laht /project/vcpkg_installed/x64-linux-dynamic/include 2>/dev/null

yum install -y curl zip unzip tar

./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install --overlay-ports=./vendor/ports --host-triplet=x64-linux-dynamic

ls -laht /project 2>/dev/null
ls -laht /project/vcpkg_installed/ 2>/dev/null
ls -laht /project/vcpkg_installed/x64-linux-dynamic 2>/dev/null
ls -laht /project/vcpkg_installed/x64-linux-dynamic/lib 2>/dev/null
ls -laht /project/vcpkg_installed/x64-linux-dynamic/include 2>/dev/null

# TODO - remove after cleanup
env

export VCPKG_ROOT="/project/vcpkg_installed/x64-linux-dynamic"
export VCPKG_ROOT="$VCPKG_ROOT"
export EIGEN_DIR="$VCPKG_ROOT/include/eigen3"
export NLOPT_DIR="$VCPKG_ROOT;$VCPKG_ROOT/lib"
export GSL_DIR="$VCPKG_ROOT;$VCPKG_ROOT/lib"
export CMAKE_BUILD_PARALLEL_LEVEL="$(nproc)"

# TODO - remove after cleanup
env
