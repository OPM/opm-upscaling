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
set(relTol 1e-4)

# Define some paths
set(RESULT_PATH ${PROJECT_BINARY_DIR}/tests/results/)
set(INPUT_DATA_PATH ${PROJECT_BINARY_DIR}/tests/input_data/)

# Create directory to store upscaling results in
file(MAKE_DIRECTORY ${RESULT_PATH})


###########################################################################
# TEST: upscale_perm on PeriodicTilted.grdecl with fixed BCs
###########################################################################

add_test(run_upscale_perm_BCp_PeriodicTilted
         ${PROJECT_BINARY_DIR}/bin/upscale_perm
	 -bc p
         -output ${RESULT_PATH}upscale_perm_BCp_PeriodicTilted.txt
         ${INPUT_DATA_PATH}grids/PeriodicTilted.grdecl)

add_test(compare_upscale_perm_BCp_PeriodicTilted
         ${PROJECT_BINARY_DIR}/bin/compare_upscaling_results
         ${INPUT_DATA_PATH}reference_solutions/upscale_perm_BCp_PeriodicTilted.txt
         ${RESULT_PATH}upscale_perm_BCp_PeriodicTilted.txt
         ${relTol}
         3 3)

set_tests_properties(compare_upscale_perm_BCp_PeriodicTilted PROPERTIES DEPENDS
                     run_upscale_perm_BCp_PeriodicTilted)
