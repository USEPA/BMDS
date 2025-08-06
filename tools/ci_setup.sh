#!/bin/bash

# print commands to stdout (x), exit if an error occurs (e)
set -xe

if [ "$RUNNER_OS" == "Linux" ]; then
    # install dependencies
    sudo apt-get update -y
    sudo apt-get install -y lcov cloc
    if [ "$generateDocx" = "true" ]; then
        sudo apt-get install -y pandoc
    fi

    # set environment variables
    export EIGEN_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-linux-dynamic/include/eigen3"
    export NLOPT_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-linux-dynamic/lib"
    export GSL_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-linux-dynamic"
    export CMAKE_BUILD_PARALLEL_LEVEL="$(nproc)"
    echo CMAKE_BUILD_PARALLEL_LEVEL="$CMAKE_BUILD_PARALLEL_LEVEL"

    # persist across Github Action steps
    echo EIGEN_DIR="$EIGEN_DIR" >> $GITHUB_ENV
    echo NLOPT_DIR="$NLOPT_DIR" >> $GITHUB_ENV
    echo GSL_DIR="$GSL_DIR" >> $GITHUB_ENV
    echo CMAKE_BUILD_PARALLEL_LEVEL="$CMAKE_BUILD_PARALLEL_LEVEL" >> $GITHUB_ENV
fi


if [ "$RUNNER_OS" == "Windows" ]; then
    # set environment variables
    export EIGEN_DIR="$GITHUB_WORKSPACE\vcpkg_installed\x64-windows\include"
    export NLOPT_DIR="$GITHUB_WORKSPACE\vcpkg_installed\x64-windows\lib"
    export GSL_DIR="$GITHUB_WORKSPACE\vcpkg_installed\x64-windows\lib"
    export CMAKE_BUILD_PARALLEL_LEVEL="$(nproc)"
    echo "CMAKE_BUILD_PARALLEL_LEVEL=$CMAKE_BUILD_PARALLEL_LEVEL"

    # persist across Github Action steps
    echo EIGEN_DIR="$EIGEN_DIR" >> $GITHUB_ENV
    echo NLOPT_DIR="$NLOPT_DIR" >> $GITHUB_ENV
    echo GSL_DIR="$GSL_DIR" >> $GITHUB_ENV
    echo CMAKE_BUILD_PARALLEL_LEVEL="$CMAKE_BUILD_PARALLEL_LEVEL" >> $GITHUB_ENV
fi


if [ "$RUNNER_OS" == "macOS" ]; then
    # set environment variables
    export EIGEN_DIR="$GITHUB_WORKSPACE/vcpkg_installed/arm64-osx-dynamic/include/eigen3"
    export NLOPT_DIR="$GITHUB_WORKSPACE/vcpkg_installed/arm64-osx-dynamic/lib"
    export GSL_DIR="$GITHUB_WORKSPACE/vcpkg_installed/arm64-osx-dynamic"
    export CMAKE_BUILD_PARALLEL_LEVEL="$(nproc)"
    echo CMAKE_BUILD_PARALLEL_LEVEL="$CMAKE_BUILD_PARALLEL_LEVEL"

    # persist across Github Action steps
    echo EIGEN_DIR="$EIGEN_DIR" >> $GITHUB_ENV
    echo NLOPT_DIR="$NLOPT_DIR" >> $GITHUB_ENV
    echo GSL_DIR="$GSL_DIR" >> $GITHUB_ENV
    echo CMAKE_BUILD_PARALLEL_LEVEL="$CMAKE_BUILD_PARALLEL_LEVEL" >> $GITHUB_ENV
fi
