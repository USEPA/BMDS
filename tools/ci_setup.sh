#!/bin/bash

if [ "$RUNNER_OS" == "Linux" ]; then

    # install dependencies
    sudo apt-get update -y
    sudo apt-get install -y automake build-essential libtool make cmake lcov cloc
    if [ "$generateDocx" = "true" ]; then
        sudo apt-get install -y pandoc
    fi

    # set export locations
    export EIGEN_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-linux-dynamic/include/eigen3"
    export NLOPT_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-linux-dynamic/lib"
    export GSL_DIR="$GITHUB_WORKSPACE/vcpkg_installed/x64-linux-dynamic"
    export CMAKE_C_COMPILER="/usr/bin/gcc-11"
    export CMAKE_CXX_COMPILER="/usr/bin/g++-11"
    export CMAKE_BUILD_PARALLEL_LEVEL="$(nproc)"
    echo "CMAKE_BUILD_PARALLEL_LEVEL=$CMAKE_BUILD_PARALLEL_LEVEL"
fi
