#!/usr/bin/env bash
set -ex

pushd . > /dev/null
opm-upscaling/travis/build-opm-upscaling.sh
cd opm-upscaling/build
ctest --output-on-failure
popd > /dev/null
