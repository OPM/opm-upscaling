# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# This file sets up five lists:
#	MAIN_SOURCE_FILES     List of compilation units which will be included in
#	                      the library. If it isn't on this list, it won't be
#	                      part of the library. Please try to keep it sorted to
#	                      maintain sanity.
#
#	TEST_SOURCE_FILES     List of programs that will be run as unit tests.
#
#	TEST_DATA_FILES       Files from the source three that should be made
#	                      available in the corresponding location in the build
#	                      tree in order to run tests there.
#
#	EXAMPLE_SOURCE_FILES  Other programs that will be compiled as part of the
#	                      build, but which is not part of the library nor is
#	                      run as tests.
#
#	PUBLIC_HEADER_FILES   List of public header files that should be
#	                      distributed together with the library. The source
#	                      files can of course include other files than these;
#	                      you should only add to this list if the *user* of
#	                      the library needs it.
#
# ATTIC_FILES           Unmaintained files. This for the projects developers
#                       only. Don't expect these files to build.

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
	)

list(APPEND MAIN_SOURCE_FILES
	opm/upscaling/RelPermUtils.cpp
	)

# originally generated with the command:
# find tests -name '*.xml' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_DATA_FILES
  tests/input_data/grids/PeriodicTilted.grdecl
  tests/input_data/grids/27cellsAniso.grdecl
  tests/input_data/grids/27cellsIso.grdecl
  tests/input_data/grids/EightCells.grdecl
  tests/input_data/grids/Hummocky.grdecl
	tests/input_data/grids/benchmark_tiny_grid.grdecl
	tests/input_data/grids/16cells.grdecl
	tests/input_data/grids/stonefile_benchmark.txt
	tests/input_data/grids/stonefileRT1.txt
	tests/input_data/grids/stonefileRT2.txt
	tests/input_data/grids/stonefileAniso.txt
	tests/input_data/grids/stonefileRT1_Thresh.txt
  tests/input_data/reference_solutions/upscale_perm_BCp_PeriodicTilted.txt 
  tests/input_data/reference_solutions/upscale_perm_BCflp_27cellsAniso.txt 
  tests/input_data/reference_solutions/upscale_perm_BCflp_27cellsIso.txt 
  tests/input_data/reference_solutions/upscale_perm_BCfl_EightCells.txt 
  tests/input_data/reference_solutions/upscale_perm_BCflp_Hummocky.txt
	tests/input_data/reference_solutions/upscale_relperm_benchmark_tiny_grid.txt
  tests/input_data/reference_solutions/upscale_relperm_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_20pnts_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_BCl_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_interpolate_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_stonefileRT1&2_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_surfaceTension_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_critRelpermThresh_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_stonefileAniso_27cellsAniso.txt
	tests/input_data/reference_solutions/upscale_relperm_gravityDailyVersion_16cells.txt
	tests/input_data/reference_solutions/upscale_relperm_gravitySBEDVersion_16cells.txt
  tests/input_data/reference_solutions/upscale_elasticity_mpc_EightCells.txt
  tests/input_data/reference_solutions/upscale_elasticity_mortar_EightCells.txt
	)

# originally generated with the command:
# find examples -name '*.c*' -a ! -name 'twophase2_test.cpp' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
	examples/cpchop.cpp
	examples/cpchop_depthtrend.cpp
	examples/cpregularize.cpp
	examples/exp_variogram.cpp
	examples/grdecldips.cpp
	examples/steadystate_test_implicit.cpp
	examples/upscale_avg.cpp
	examples/upscale_cap.cpp
	examples/upscale_cond.cpp
	examples/upscale_perm.cpp
	examples/upscale_relperm.cpp
	examples/upscale_relpermvisc.cpp
	examples/upscale_singlephase.cpp
	examples/upscale_steadystate_implicit.cpp
	tests/compare_upscaling_results.cpp
	)

list (APPEND ADDITIONAL_SOURCE_FILES
  benchmarks/upscale_relperm_benchmark.cpp
	   )

# originally generated with the command:
# find attic -name '*.c*' -printf '\t%p\n' | sort
list (APPEND ATTIC_FILES
	attic/aniso_implicit_steadystate_test.cpp
	attic/aniso_steadystate_test.cpp
	attic/implicit_steadystate_test.cpp
	attic/steadystate_test_explicit.cpp
	)

# programs listed here will not only be compiled, but also marked for
# installation
list (APPEND PROGRAM_SOURCE_FILES
	examples/cpchop.cpp
	examples/cpchop_depthtrend.cpp
	examples/cpregularize.cpp
	examples/exp_variogram.cpp
	examples/grdecldips.cpp
	examples/steadystate_test_implicit.cpp
	examples/upscale_avg.cpp
	examples/upscale_cap.cpp
	examples/upscale_cond.cpp
	examples/upscale_perm.cpp
	examples/upscale_relperm.cpp
	examples/upscale_relpermvisc.cpp
	examples/upscale_singlephase.cpp
	examples/upscale_steadystate_implicit.cpp
	)

# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
	opm/upscaling/ParserAdditions.hpp
	opm/upscaling/SinglePhaseUpscaler.hpp
	opm/upscaling/SteadyStateUpscaler.hpp
	opm/upscaling/SteadyStateUpscaler_impl.hpp
	opm/upscaling/SteadyStateUpscalerImplicit.hpp
	opm/upscaling/SteadyStateUpscalerImplicit_impl.hpp
	opm/upscaling/SteadyStateUpscalerManager.hpp
	opm/upscaling/SteadyStateUpscalerManagerImplicit.hpp
	opm/upscaling/UpscalerBase.hpp
	opm/upscaling/UpscalerBase_impl.hpp
	opm/upscaling/UpscalingTraits.hpp
	)
