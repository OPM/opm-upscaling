#!/bin/bash

declare -a upstreams
upstreams=(opm-common
           opm-grid)

declare -A upstreamRev
upstreamRev[opm-common]=master
upstreamRev[opm-grid]=master

# Currently no downstreams
declare -a downstreams
declare -A downstreamRev

# Fetch opm-common before loading its shared Jenkins helpers.
source "$WORKSPACE/jenkins/checkout-opm-common.sh"

source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

parseRevisions
printHeader opm-upscaling

clone_repositories opm-upscaling

build_module_full opm-upscaling
