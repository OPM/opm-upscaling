#!/bin/bash

declare -a upstreams
upstreams=(opm-common
           opm-grid)

declare -A upstreamRev
upstreamRev[opm-common]=master
upstreamRev[opm-grid]=master

# Fetch opm-common before loading its shared Jenkins helpers.
source "$WORKSPACE/jenkins/checkout-opm-common.sh"

source ${TOOLCHAIN_DIR}/build-configurations-sca.sh
source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

$WORKSPACE/jenkins/build.sh

run_static_analysis opm-upscaling
