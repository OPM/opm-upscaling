#!/bin/bash

# Clear output directories
rm -rf inversion_output
mkdir inversion_output

# Run simulator
../../../../../examples/blackoil_sim_test inversiontest.param > inversion_output/stdout_dump 2>&1

# For Atgeirr only. ParaView on mac crashes on denormal floating point numbers, below fixes.
# fsubst0 Float32 Float64 gravity_output/*.vtu > /dev/null
#rm gravity_output/*~~

# Pruning output
rm -f inversion_output/blackoil-output*[1-9].vtu
rm -f inversion_output/blackoil-output*[1-9].dat

# If this diff only gives timing differences, we are happy.
diff inversion_reference_output inversion_output
