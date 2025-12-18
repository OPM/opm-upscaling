#!/bin/bash

declare -a upstreams
upstreams=(opm-common
           opm-grid)

declare -A upstreamRev
upstreamRev[opm-common]=master
upstreamRev[opm-grid]=master

if grep -q "opm-common=" <<< $ghprbCommentBody
then
  upstreamRev[opm-common]=pull/`echo $ghprbCommentBody | sed -r 's/.*opm-common=([0-9]+).*/\1/g'`/merge
fi

# Currently no downstreams
declare -a downstreams
declare -A downstreamRev

# Clone opm-common
if ! test -d $WORKSPACE/deps/opm-common
then
  mkdir -p $WORKSPACE/deps/opm-common
  pushd $WORKSPACE/deps/opm-common
  git init .
  git remote add origin https://github.com/OPM/opm-common
  git fetch --depth 1 origin ${upstreamRev[opm-common]}:branch_to_build
  test $? -eq 0 || exit 1
  git checkout branch_to_build
  popd
fi

source $WORKSPACE/deps/opm-common/jenkins/build-opm-module.sh

parseRevisions
printHeader opm-upscaling

clone_repositories opm-upscaling

build_module_full opm-upscaling
