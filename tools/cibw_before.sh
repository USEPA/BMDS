#!/bin/bash

set -eo pipefail

if [ "${RUNNER_OS}" = "Linux" ]; then
    yum install -y zip
fi

# Dynamically locate vcpkg inside or outside the mapped repo root
if [ -d "./vcpkg" ]; then
    VCPKG_PATH="./vcpkg"
elif [ -d "../vcpkg" ]; then
    VCPKG_PATH="../vcpkg"
else
    echo "ERROR: vcpkg directory not found!"
    exit 1
fi

echo "Found vcpkg at: $VCPKG_PATH"
echo "Host Triplet to build: $VCPKG_HOST_TRIPLET"

./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install --host-triplet=$VCPKG_HOST_TRIPLET --x-manifest-root=.

pip install pybind11==3.0.0 --target=./pybind11
