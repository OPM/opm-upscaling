#!/bin/bash
set -e

INPUT_DATA_PATH="$1"
RESULT_PATH="$2"
BINPATH="$3"
TEST_NAME="$4"
ABS_TOL="$5"
REL_TOL="$6"
EXE_NAME="$7"
shift 7
TEST_ARGS="$@"

rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}

${BINPATH}/${EXE_NAME} ${TEST_ARGS}
${BINPATH}/compareUpscaling ${INPUT_DATA_PATH}/reference_solutions/${TEST_NAME}.txt ${RESULT_PATH}/${TEST_NAME}.txt ${ABS_TOL} ${REL_TOL}

rm -Rf ${RESULT_PATH}
