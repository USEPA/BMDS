#!/bin/bash

set -xe

cd $GITHUB_WORKSPACE

brew install gsl --quiet
brew install --formula ./tools/mac-brew-nlopt-271.rb --quiet
brew install eigen --quiet
