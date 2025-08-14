#!/bin/bash

set -eo pipefail

if [ "${RUNNER_OS}" = "Linux" ]; then
    yum install -y zip
fi

echo $VCPKG_HOST_TRIPLET
./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install --host-triplet=$VCPKG_HOST_TRIPLET

pip install pybind11==3.0.0 --target=./pybind11
