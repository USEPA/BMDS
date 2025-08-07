#!/bin/bash

# exit if an error occurs
set -xe

env

pwd

ls -laht /project | true
ls -laht /project/vcpkg_installed/ | true
ls -laht /project/vcpkg_installed/x64-linux-dynamic | true
ls -laht /project/vcpkg_installed/x64-linux-dynamic/lib | true
ls -laht /project/vcpkg_installed/x64-linux-dynamic/include | true

yum update -y
yum install -y zip

./vcpkg/bootstrap-vcpkg.sh
./vcpkg/vcpkg install --overlay-ports=./vendor/ports --host-triplet=x64-linux-dynamic
