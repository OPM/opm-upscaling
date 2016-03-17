# Add tests to check if the upscaling binaries calculates correctly
#
# Tests are added in two steps (see already added tests below for examples):
#
#    1) Add test that runs the binary and output results to a file:
#       add_test(<testname> <command>)
#       <command> refers to the command used to run the binary with input variables in a terminal
#
#    2) Add test that compares the output from the previous test:
#       add_test(<testname> ${PROJECT_BINARY_DIR}/bin/compare_upscaling_results <path_to_refSoln> <path_to_newSoln>
#                ${tol} <number_of_result_rows> <number_of_result_cols>)
#       This test should depend on the first test, so include:
#       set_tests_properties(<test1> PROPERTIES DEPENDS <test2>)
#
# Some naming conventions:
#    The first test should be named:  run_<binary_name>_<options>_<model_name>
#    The second test:                 compare_<binary_name>_<options>_<model_name>
#
# Test models and reference solutions are available in ${PROJECT_BINARY_DIR}/tests/input_data.
# New test data can be made available in the build tree by including them in CMakeLists_files.cmake,
# under 'APPEND TEST_SOURCE_FILES'.


# Set absolute tolerance to be used for testing
set(tol 1e-2)

# Define some paths
set(RESULT_PATH ${PROJECT_BINARY_DIR}/tests/results)
set(INPUT_DATA_PATH ${PROJECT_BINARY_DIR}/tests/input_data)

# Create directory to store upscaling results in
file(MAKE_DIRECTORY ${RESULT_PATH})


###########################################################################
# TEST: upscale_perm 
###########################################################################

# Define macro that performs the two steps mentioned above for upscale_perm
# Input: 
#   - gridname: basename (no extension) of grid model
#   - bcs: Boundary condition type (f, l or p, or combinations of these)
#   - rows: Number of rows in result file that is to be compared
# This macro assumes that ${gridname}.grdecl is found in directory ${INPUT_DATA_PATH}grids/
# and that upscale_perm_BC${bcs}_${gridname}.txt is found in ${INPUT_DATA_PATH}reference_solutions
macro (add_test_upscale_perm gridname bcs rows)
  # Add test that runs upscale_perm and outputs the results to file
  opm_add_test(upscale_perm_BC${bcs}_${gridname} NO_COMPILE
               EXE_NAME upscale_perm
               DRIVER_ARGS ${INPUT_DATA_PATH} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           upscale_perm_BC${bcs}_${gridname}
                           ${tol} ${rows} 3
               TEST_ARGS -bc ${bcs}
                         -output ${RESULT_PATH}/upscale_perm_BC${bcs}_${gridname}.txt
                         ${INPUT_DATA_PATH}/grids/${gridname}.grdecl)
endmacro (add_test_upscale_perm gridname bcs)

###########################################################################
# TEST: upscale_relperm 
###########################################################################

# Define macro that performs the two steps mentioned above for upscale_relperm
# Input: 
#   - gridname: basename (no extension) of grid model
#   - rows: Number of rows in result file that is to be compared
# This macro assumes that ${gridname}.grdecl is found in directory ${INPUT_DATA_PATH}grids/
# and that upscale_perm_BC${bcs}_${gridname}.txt is found in ${INPUT_DATA_PATH}reference_solutions
macro (add_test_upscale_relperm testname gridname stonefiles rows cols)
  # Add test that runs upscale_perm and outputs the results to file
  set(test_args ${ARGN}
                -output ${RESULT_PATH}/upscale_relperm_${testname}.txt
                ${INPUT_DATA_PATH}/grids/${gridname}.grdecl)
  foreach(stonefile ${stonefiles})
    list(APPEND test_args ${INPUT_DATA_PATH}/grids/${stonefile})
  endforeach()
  opm_add_test(upscale_relperm_${testname} NO_COMPILE
               EXE_NAME upscale_relperm
               DRIVER_ARGS ${INPUT_DATA_PATH} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           upscale_relperm_${testname}
                           0.02 ${rows} ${cols}
               TEST_ARGS ${test_args})
endmacro ()

###########################################################################
# TEST: upscale_elasticity
###########################################################################

# Define macro that performs the two steps mentioned above for upscale_elasticity
# Input: 
#   - gridname: basename (no extension) of grid model
#   - method: method to apply
# This macro assumes that ${gridname}.grdecl is found in directory ${INPUT_DATA_PATH}grids/
# and that upscale_elasticity_${method}_${gridname}.txt is found in ${INPUT_DATA_PATH}reference_solutions
macro (add_test_upscale_elasticity gridname method)
  # Add test that runs upscale_perm and outputs the results to file
  opm_add_test(upscale_elasticity_${method}_${gridname} NO_COMPILE
               EXE_NAME upscale_elasticity
               DRIVER_ARGS ${INPUT_DATA_PATH} ${RESULT_PATH}
                           ${CMAKE_BINARY_DIR}/bin
                           upscale_elasticity_${method}_${gridname}
                           ${tol} 6 6
               TEST_ARGS output=${RESULT_PATH}/upscale_elasticity_${method}_${gridname}.txt
                         gridfilename=${INPUT_DATA_PATH}/grids/${gridname}.grdecl
                         method=${method})
endmacro (add_test_upscale_elasticity gridname method rows)

# Make sure that we build the helper executable before running tests
# (the "tests" target is setup in OpmLibMain.cmake)
if(NOT TARGET test-suite)
  add_custom_target(test-suite)
endif()
add_dependencies (test-suite datafiles upscale_perm upscale_relperm)
add_dependencies (test-suite compare_upscaling_results)
opm_set_test_driver(${PROJECT_SOURCE_DIR}/tests/runtest.sh "")

# Add tests for different models
add_test_upscale_perm(PeriodicTilted p 3)
add_test_upscale_perm(27cellsAniso flp 9)
add_test_upscale_perm(27cellsIso flp 9)
add_test_upscale_perm(EightCells fl 6)
add_test_upscale_perm(Hummocky flp 9)

# Add tests for different models
add_test_upscale_relperm(BCf_pts20_surfTens11_stonefile_benchmark_stonefile_benchmark_benchmark_tiny_grid
                         benchmark_tiny_grid stonefile_benchmark.txt 20 8
                         -bc f -points 20 -relPermCurve 2 -upscaleBothPhases true
                         -jFunctionCurve 3 -surfaceTension 11 -gravity 0.0
                         -waterDensity 1.0 -oilDensity 0.6 -interpolate 0
                         -maxpoints 1000 -outputprecision 20 -maxPermContrast 1e7
                         -minPerm 1e-12 -maxPerm 100000 -minPoro 0.0001
                         -saturationThreshold 0.0001 -linsolver_tolerance 1e-12
                         -linsolver_verbosity 0 -linsolver_type 3 -fluids ow
                         -krowxswirr -1 -krowyswirr -1 -krowzswirr -1
                         -doEclipseCheck true -critRelpermThresh 1e-6)
add_test_upscale_relperm(BCf_pts30_surfTens11_stone1_stone1_EightCells
                         EightCells stone1.txt 30 8)
add_test_upscale_relperm(BCf_pts20_surfTens11_stone1_stone1_EightCells
                         EightCells stone1.txt 30 8 -points 20)
add_test_upscale_relperm(BCl_pts30_surfTens11_stone1_stone1_EightCells
                         EightCells stone1.txt 30 20 -bc l)
add_test_upscale_relperm(BCf_pts30_surfTens45_stone1_stone1_EightCells
                         EightCells stone1.txt 30 8 -surfaceTension 45)
add_test_upscale_relperm(BCf_pts30_surfTens11_stone1_stone2_EightCells
                         EightCells "stone1.txt;stone2.txt" 30 8)
add_test_upscale_relperm(BCf_pts30_surfTens11_stoneAniso_stoneAniso_27cellsAniso
                         27cellsAniso stoneAniso.txt 30 8)

if((DUNE_ISTL_VERSION_MAJOR GREATER 2) OR
   (DUNE_ISTL_VERSION_MAJOR EQUAL 2 AND DUNE_ISTL_VERSION_MINOR GREATER 2))
  add_dependencies (test-suite upscale_elasticity)
  add_test_upscale_elasticity(EightCells mpc)
  add_test_upscale_elasticity(EightCells mortar)
endif()

add_custom_target(clear_test_data ${CMAKE_COMMAND} -E remove_directory ${RESULT_PATH}
                  COMMAND ${CMAKE_COMMAND} -E make_directory ${RESULT_PATH}
                  COMMENT "Removing old test results")

add_dependencies(check clear_test_data)
