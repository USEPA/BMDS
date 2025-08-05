#!/bin/bash

# print commands to stdout (x), exit if an error occurs (e)
set -xe

if [ "$RUNNER_OS" == "Linux" ]; then

    # install dependencies
    sudo apt-get update -y
    sudo apt-get install -y automake build-essential libtool make cmake lcov cloc
    if [ "$generateDocx" = "true" ]; then
        sudo apt-get install -y pandoc
    fi

    # set environment variables
    export EIGEN_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-linux-dynamic/include/eigen3"
    export NLOPT_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-linux-dynamic/lib"
    export GSL_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-linux-dynamic"
    export CMAKE_C_COMPILER="/usr/bin/gcc-11"
    export CMAKE_CXX_COMPILER="/usr/bin/g++-11"
    export CMAKE_BUILD_PARALLEL_LEVEL="$(nproc)"
    echo "CMAKE_BUILD_PARALLEL_LEVEL=$CMAKE_BUILD_PARALLEL_LEVEL"

    # persist across Github Action steps
    echo "EIGEN_DIR=$EIGEN_DIR" >> $GITHUB_ENV
    echo "NLOPT_DIR=$NLOPT_DIR" >> $GITHUB_ENV
    echo "GSL_DIR=$GSL_DIR" >> $GITHUB_ENV
    echo "CMAKE_C_COMPILER=$CMAKE_C_COMPILER" >> $GITHUB_ENV
    echo "CMAKE_CXX_COMPILER=$CMAKE_CXX_COMPILER" >> $GITHUB_ENV
    echo "CMAKE_BUILD_PARALLEL_LEVEL=$CMAKE_BUILD_PARALLEL_LEVEL" >> $GITHUB_ENV
fi


if [ "$RUNNER_OS" == "Windows" ]; then
    # set environment variables
    export EIGEN_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-windows-static/include"
    export NLOPT_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-windows-static"
    export GSL_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-windows-static"
    export CMAKE_BUILD_PARALLEL_LEVEL="$(nproc)"
    echo "CMAKE_BUILD_PARALLEL_LEVEL=$CMAKE_BUILD_PARALLEL_LEVEL"

    # persist across Github Action steps
    echo "EIGEN_DIR=$EIGEN_DIR" >> $GITHUB_ENV
    echo "NLOPT_DIR=$NLOPT_DIR" >> $GITHUB_ENV
    echo "GSL_DIR=$GSL_DIR" >> $GITHUB_ENV
    echo "CMAKE_BUILD_PARALLEL_LEVEL=$CMAKE_BUILD_PARALLEL_LEVEL" >> $GITHUB_ENV
fi


if [ "$RUNNER_OS" == "macOS" ]; then
    echo "Not implemented (yet)"
    exit 1
fi
