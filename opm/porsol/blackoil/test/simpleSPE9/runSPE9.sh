#!/bin/bash

# Clear output directories
rm -rf spe9_output
mkdir spe9_output

# Run simulator
../../../../../examples/blackoil_sim_test spe9run6500days.xml > spe9_output/stdout_dump 2>&1

# Pruning output
rm -f spe9_output/blackoil-output*[1-9].vtu
rm -f spe9_output/blackoil-output*[1-9].dat
rm -f spe9_output/blackoil-output*[1-9]0.vtu
rm -f spe9_output/blackoil-output*[1-9]0.dat


# If this diff only gives timing differences, we are happy.
diff spe9_reference_output spe9_output
