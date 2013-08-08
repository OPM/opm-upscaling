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
# find opm -name '*.c*' -printf '\t%p\n' | sort
list (APPEND MAIN_SOURCE_FILES
	opm/elasticity/boundarygrid.cpp
	opm/elasticity/material.cpp
	opm/elasticity/materials.cpp
	opm/elasticity/matrixops.cpp
	opm/elasticity/mpc.cpp
	)

# originally generated with the command:
# find tests -name '*.cpp' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_SOURCE_FILES
	)

# originally generated with the command:
# find tests -name '*.xml' -a ! -wholename '*/not-unit/*' -printf '\t%p\n' | sort
list (APPEND TEST_DATA_FILES
  tests/input_data/grids/PeriodicTilted.grdecl
  tests/input_data/grids/27cellsAniso.grdecl
  tests/input_data/grids/27cellsIso.grdecl
  tests/input_data/grids/EightCells.grdecl
  tests/input_data/grids/Hummocky.grdecl
  tests/input_data/reference_solutions/upscale_perm_BCp_PeriodicTilted.txt 
  tests/input_data/reference_solutions/upscale_perm_BCflp_27cellsAniso.txt 
  tests/input_data/reference_solutions/upscale_perm_BCflp_27cellsIso.txt 
  tests/input_data/reference_solutions/upscale_perm_BCfl_EightCells.txt 
  tests/input_data/reference_solutions/upscale_perm_BCflp_Hummocky.txt 
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
	examples/upscale_elasticity.cpp
	examples/upscale_perm.cpp
	examples/upscale_relperm.cpp
	examples/upscale_relpermvisc.cpp
	examples/upscale_singlephase.cpp
	examples/upscale_steadystate_implicit.cpp
	tests/compare_upscaling_results.cpp
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
	examples/upscale_elasticity.cpp
	examples/upscale_perm.cpp
	examples/upscale_relperm.cpp
	examples/upscale_relpermvisc.cpp
	examples/upscale_singlephase.cpp
	examples/upscale_steadystate_implicit.cpp
	)

# originally generated with the command:
# find opm -name '*.h*' -a ! -name '*-pch.hpp' -printf '\t%p\n' | sort
list (APPEND PUBLIC_HEADER_FILES
	opm/elasticity/applier.hpp
	opm/elasticity/asmhandler.hpp
	opm/elasticity/asmhandler_impl.hpp
	opm/elasticity/boundarygrid.hh
	opm/elasticity/elasticity.hpp
	opm/elasticity/elasticity_impl.hpp
	opm/elasticity/elasticity_upscale.hpp
	opm/elasticity/elasticity_upscale_impl.hpp
	opm/elasticity/logutils.hpp
	opm/elasticity/material.hh
	opm/elasticity/materials.hh
	opm/elasticity/matrixops.hpp
	opm/elasticity/mortar_evaluator.hpp
	opm/elasticity/mortar_schur.hpp
	opm/elasticity/mortar_schur_precond.hpp
	opm/elasticity/mortar_utils.hpp
	opm/elasticity/mpc.hh
	opm/elasticity/seqlu.hpp
	opm/elasticity/shapefunctions.hpp
	opm/elasticity/uzawa_solver.hpp
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
