#!/bin/bash

set -xe

cd $GITHUB_WORKSPACE

brew install --formula ./tools/mac/gsl.rb --quiet
brew install --formula ./tools/mac/nlopt.rb --quiet
brew install --formula ./tools/mac/eigen.rb --quiet
