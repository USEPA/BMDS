#!/bin/bash

# exit if an error occurs
set -e

env

if [ "$RUNNER_OS" == "Linux" ]; then
./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install --overlay-ports="./vendor/ports" --host-triplet="x64-linux-dynamic"
fi


if [ "$RUNNER_OS" == "Windows" ]; then

fi


if [ "$RUNNER_OS" == "macOS" ]; then

fi
