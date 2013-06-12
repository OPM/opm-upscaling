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
#                ${relTol} <number_of_result_rows> <number_of_result_cols>)
#       This test should depend on the first test, so include:
#       set_tests_properties(<test1> PROPERTIES DEPENDS <test2>)
#
# Some naming conventions:
#    The first test should be named:  run_<binary_name>_<options>_<model_name>
#    The second test:                 compare_<binary_name>_<options>_<model_name>
#
# Test models and reference solutions are available in ${PROJECT_BINARY_DIR}/tests/input_data.
# New test data can be made available in the build three by including them in CMakeLists_files.cmake,
# under 'APPEND TEST_SOURCE_FILES'.


# Set relative tolerance to be used for testing
set(relTol 1e-2)

# Define some paths
set(RESULT_PATH ${PROJECT_BINARY_DIR}/tests/results/)
set(INPUT_DATA_PATH ${PROJECT_BINARY_DIR}/tests/input_data/)

# Create directory to store upscaling results in
file(MAKE_DIRECTORY ${RESULT_PATH})


###########################################################################
# TEST: upscale_perm 
###########################################################################

# Define macro that performs the two steps mentioned above for upscale_perm
# Input: 
#   - gridname: basename (no extension) of grid model
#   - bcs: Boundary condition type (f, l or p or combinations)
#   - rows: Number of rows in result file that is to be compared
# This macro assumes that ${gridname}.grdecl is found in directory ${INPUT_DATA_PATH}grids/
# and that upscale_perm_BC${bcs}_${gridname}.txt is found in ${INPUT_DATA_PATH}reference_solutions
macro (add_test_upscale_perm gridname bcs rows)
  # Add test that runs upscale_perm and outputs the results to file
  add_test(run_upscale_perm_BC${bcs}_${gridname}
	   ${PROJECT_BINARY_DIR}/bin/upscale_perm
	   -bc ${bcs}
	   -output ${RESULT_PATH}upscale_perm_BC${bcs}_${gridname}.txt
	   ${INPUT_DATA_PATH}grids/${gridname}.grdecl)
  # Add test that compare the results from the previous test with a reference solution
  add_test(compare_upscale_perm_BC${bcs}_${gridname}
           ${PROJECT_BINARY_DIR}/bin/compare_upscaling_results
           ${INPUT_DATA_PATH}reference_solutions/upscale_perm_BC${bcs}_${gridname}.txt
           ${RESULT_PATH}upscale_perm_BC${bcs}_${gridname}.txt
           ${relTol}
           ${rows} 3)
  # Set dependency of the two tests
  set_tests_properties(compare_upscale_perm_BC${bcs}_${gridname} PROPERTIES DEPENDS
                       run_upscale_perm_BC${bcs}_${gridname})
endmacro (add_test_upscale_perm gridname bcs)

# Add tests for different models
add_test_upscale_perm(PeriodicTilted p 3)
add_test_upscale_perm(27cellsAniso flp 9)
add_test_upscale_perm(27cellsIso flp 9)
add_test_upscale_perm(EightCells fl 6)
add_test_upscale_perm(hummockyChopped flp 9)
