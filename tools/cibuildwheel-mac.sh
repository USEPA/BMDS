#!/bin/bash

set -xe

brew install gsl --quiet
brew install nlopt --quiet
brew install eigen --quiet

cd $GITHUB_WORKSPACE
