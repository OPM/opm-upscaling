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
	tests/common/boundaryconditions_test.cpp
	tests/common/matrix_test.cpp
	tests/common/test_gravitypressure.cpp
	)

list(APPEND MAIN_SOURCE_FILES
	opm/porsol/blackoil/fluid/BlackoilPVT.cpp
	opm/porsol/blackoil/fluid/MiscibilityDead.cpp
	opm/porsol/blackoil/fluid/MiscibilityLiveGas.cpp
	opm/porsol/blackoil/fluid/MiscibilityLiveOil.cpp
	opm/porsol/blackoil/fluid/MiscibilityProps.cpp
	opm/porsol/common/blas_lapack.cpp
	opm/porsol/common/BoundaryPeriodicity.cpp
	opm/porsol/common/ImplicitTransportDefs.cpp
	opm/porsol/common/setupGridAndProps.cpp
	opm/porsol/euler/ImplicitCapillarity.cpp
	opm/upscaling/ParserAdditions.cpp
	opm/upscaling/RelPermUtils.cpp
        opm/upscaling/initCPGrid.cpp
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
	tests/input_data/grids/stone1.txt
	tests/input_data/grids/stone2.txt
	tests/input_data/grids/stoneAniso.txt
	tests/input_data/grids/stonefile_benchmark.txt
  tests/input_data/reference_solutions/upscale_perm_BCp_PeriodicTilted.txt 
  tests/input_data/reference_solutions/upscale_perm_BCflp_27cellsAniso.txt 
  tests/input_data/reference_solutions/upscale_perm_BCflp_27cellsIso.txt 
  tests/input_data/reference_solutions/upscale_perm_BCfl_EightCells.txt 
  tests/input_data/reference_solutions/upscale_perm_BCflp_Hummocky.txt
	tests/input_data/reference_solutions/upscale_relperm_BCf_pts20_surfTens11_stonefile_benchmark_stonefile_benchmark_benchmark_tiny_grid.txt
	tests/input_data/reference_solutions/upscale_relperm_BCf_pts30_surfTens11_stone1_stone1_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_BCf_pts30_surfTens11_stone1_stone2_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_BCf_pts20_surfTens11_stone1_stone1_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_BCl_pts30_surfTens11_stone1_stone1_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_BCf_pts30_surfTens45_stone1_stone1_EightCells.txt
	tests/input_data/reference_solutions/upscale_relperm_BCf_pts30_surfTens11_stoneAniso_stoneAniso_27cellsAniso.txt
  tests/input_data/reference_solutions/upscale_elasticity_mpc_EightCells.txt
  tests/input_data/reference_solutions/upscale_elasticity_mortar_EightCells.txt
	)

# originally generated with the command:
# find examples -name '*.c*' -a ! -name 'twophase2_test.cpp' -printf '\t%p\n' | sort
list (APPEND EXAMPLE_SOURCE_FILES
	examples/aniso_implicitcap_test.cpp
	examples/aniso_simulator_test.cpp
	examples/co2_blackoil_pvt.cpp
	examples/cpchop.cpp
	examples/cpchop_depthtrend.cpp
	examples/cpregularize.cpp
	examples/exp_variogram.cpp
	examples/grdecldips.cpp
	examples/implicitcap_test.cpp
	examples/known_answer_test.cpp
	examples/mimetic_aniso_solver_test.cpp
	examples/mimetic_periodic_test.cpp
	examples/mimetic_solver_test.cpp
	examples/sim_blackoil_impes.cpp
	examples/sim_co2_impes.cpp
	examples/sim_steadystate_explicit.cpp
	examples/sim_steadystate_implicit.cpp
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
	examples/sim_blackoil_impes.cpp
	examples/sim_co2_impes.cpp
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
	opm/porsol/blackoil/BlackoilFluid.hpp
	opm/porsol/blackoil/BlackoilInitialization.hpp
	opm/porsol/blackoil/BlackoilSimulator.hpp
	opm/porsol/blackoil/BlackoilWells.hpp
	opm/porsol/blackoil/co2fluid/benchmark3co2tables.hh
	opm/porsol/blackoil/co2fluid/BlackoilCo2PVT.hpp
	opm/porsol/blackoil/ComponentTransport.hpp
	opm/porsol/blackoil/fluid/BlackoilComponent.hpp
	opm/porsol/blackoil/fluid/BlackoilDefs.hpp
	opm/porsol/blackoil/fluid/BlackoilPVT.hpp
	opm/porsol/blackoil/fluid/FluidMatrixInteractionBlackoil.hpp
	opm/porsol/blackoil/fluid/FluidStateBlackoil.hpp
	opm/porsol/blackoil/fluid/MiscibilityDead.hpp
	opm/porsol/blackoil/fluid/MiscibilityLiveGas.hpp
	opm/porsol/blackoil/fluid/MiscibilityLiveOil.hpp
	opm/porsol/blackoil/fluid/MiscibilityProps.hpp
	opm/porsol/blackoil/fluid/MiscibilityWater.hpp
	opm/porsol/common/BCRSMatrixBlockAssembler.hpp
	opm/porsol/common/blas_lapack.hpp
	opm/porsol/common/BoundaryConditions.hpp
	opm/porsol/common/BoundaryPeriodicity.hpp
	opm/porsol/common/fortran.hpp
	opm/porsol/common/GridInterfaceEuler.hpp
	opm/porsol/common/ImplicitTransportDefs.hpp
	opm/porsol/common/Matrix.hpp
	opm/porsol/common/MatrixInverse.hpp
	opm/porsol/common/PeriodicHelpers.hpp
	opm/porsol/common/ReservoirPropertyCapillaryAnisotropicRelperm.hpp
	opm/porsol/common/ReservoirPropertyCapillaryAnisotropicRelperm_impl.hpp
	opm/porsol/common/ReservoirPropertyCapillary.hpp
	opm/porsol/common/ReservoirPropertyCapillary_impl.hpp
	opm/porsol/common/ReservoirPropertyCommon.hpp
	opm/porsol/common/ReservoirPropertyCommon_impl.hpp
	opm/porsol/common/ReservoirPropertyFixedMobility.hpp
	opm/porsol/common/ReservoirPropertyTracerFluid.hpp
	opm/porsol/common/RockAnisotropicRelperm.hpp
	opm/porsol/common/Rock.hpp
	opm/porsol/common/Rock_impl.hpp
	opm/porsol/common/RockJfunc.hpp
	opm/porsol/common/setupBoundaryConditions.hpp
	opm/porsol/common/setupGridAndProps.hpp
	opm/porsol/common/SimulatorBase.hpp
	opm/porsol/common/SimulatorTraits.hpp
	opm/porsol/common/SimulatorUtilities.hpp
	opm/porsol/common/Wells.hpp
	opm/porsol/euler/CflCalculator.hpp
	opm/porsol/euler/EulerUpstream.hpp
	opm/porsol/euler/EulerUpstream_impl.hpp
	opm/porsol/euler/EulerUpstreamImplicit.hpp
	opm/porsol/euler/EulerUpstreamImplicit_impl.hpp
	opm/porsol/euler/EulerUpstreamResidual.hpp
	opm/porsol/euler/EulerUpstreamResidual_impl.hpp
	opm/porsol/euler/ImplicitCapillarity.hpp
	opm/porsol/euler/ImplicitCapillarity_impl.hpp
	opm/porsol/euler/MatchSaturatedVolumeFunctor.hpp
	opm/porsol/mimetic/IncompFlowSolverHybrid.hpp
	opm/porsol/mimetic/MimeticIPAnisoRelpermEvaluator.hpp
	opm/porsol/mimetic/MimeticIPEvaluator.hpp
	opm/porsol/mimetic/TpfaCompressibleAssembler.hpp
	opm/porsol/mimetic/TpfaCompressible.hpp
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
        opm/upscaling/initCPGrid.hpp
	)
