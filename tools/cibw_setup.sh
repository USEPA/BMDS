#!/bin/bash

# exit if an error occurs
set -e

env

./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install --overlay-ports="./vendor/ports" --host-triplet="x64-linux-dynamic"
