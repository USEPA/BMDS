#!/bin/bash

# exit if an error occurs
set -e

echo "Current OS: $RUNNER_OS"
echo "$GITHUB_WORKSPACE"
ls $GITHUB_WORKSPACE

ls $GITHUB_WORKSPACE/vcpkg || true
ls $GITHUB_WORKSPACE/vcpkg_installed || true

if [ "$RUNNER_OS" == "Linux" ]; then
$GITHUB_WORKSPACE/vcpkg/bootstrap-vcpkg.sh
$GITHUB_WORKSPACE/vcpkg/vcpkg install --overlay-ports="./vendor/ports" --host-triplet="x64-linux-dynamic"
fi


if [ "$RUNNER_OS" == "Windows" ]; then

fi


if [ "$RUNNER_OS" == "macOS" ]; then

fi
