#!/bin/bash

source `dirname $0`/build-opm-upscaling.sh

ERT_REVISION=master
OPM_COMMON_REVISION=master
OPM_PARSER_REVISION=master
OPM_MATERIAL_REVISION=master
OPM_CORE_REVISION=master
DUNE_CORNERPOINT_REVISION=master
OPM_OUTPUT_REVISION=master

build_opm_upscaling
test $? -eq 0 || exit 1

cp serial/build-opm-upscaling/testoutput.xml .
