# Create directory to store upscaling results
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/tests/results)

add_test(run_upscale_perm 
	 ${PROJECT_BINARY_DIR}/bin/upscale_perm 
	 -output ${PROJECT_BINARY_DIR}/tests/results/upscale_perm_PeriodicTilted.txt
         ${PROJECT_BINARY_DIR}/tests/input_data/grids/PeriodicTilted.grdecl)