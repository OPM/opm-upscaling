# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

# set up project and specify the minimum cmake version
project("opm-upscaling" C CXX)

cmake_minimum_required(VERSION 3.0)
option(INSTALL_BENCHMARKS "Install benchmark applications?" OFF)

# Set CMP0053 (how to handle escape sequences in strings) to the new
# behavior to avoid a pretty annoying cmake warning if a library is
# defined in the toplevel CMakeLists.txt. This should probably be
# considered to be a bug in the dune build system. Note that the old
# behaviour will most likely also work fine, but the result of setting
# this policy to NEW is most likely what's intended.
if (POLICY CMP0053)
  cmake_policy(SET CMP0053 NEW)
endif()

# find the build system (i.e., dune-common) and set cmake's module path
find_package(dune-common REQUIRED)
list(APPEND CMAKE_MODULE_PATH ${dune-common_MODULE_PATH}
  "${PROJECT_SOURCE_DIR}/cmake/modules")

# include the dune macros
include(DuneMacros)

# start a dune project with information from dune.module
dune_project()

# include the OPM cmake macros
include(OpmMacros)

# find the packages needed to compile the module
find_package(Boost COMPONENTS unit_test_framework iostreams REQUIRED)

# recursively mark all header files beneath the "opm" directory for
# installation.
opm_recusive_export_all_headers("opm")

# the cmake modules get a special treatment
opm_export_cmake_modules()

# we want all features detected by the build system to be enabled,
# thank you!
dune_enable_all_packages()

# add libopmupscaling.a and specify that it encompasses all source files
# which are located somewhere beneath the "opm" directory. Note that
# if the tests and the source files of the library were intermingled,
# this would be not so easy, i.e. it would require a separate
# CMakeLists.txt file in each subdirectory and a call to
# add_subdirectory() for each of them.
opm_recursive_add_library("opmupscaling" "opm")

# setup the tests
opm_add_test(boundaryconditions_test
  SOURCES tests/common/boundaryconditions_test.cpp
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(matrix_test
  SOURCES tests/common/matrix_test.cpp
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_test(test_gravitypressure
  SOURCES tests/common/test_gravitypressure.cpp
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")

opm_add_application(aniso_implicitcap_test SOURCES examples/aniso_implicitcap_test.cpp)
opm_add_application(aniso_simulator_test SOURCES examples/aniso_simulator_test.cpp)
opm_add_application(co2_blackoil_pvt SOURCES examples/co2_blackoil_pvt.cpp)
opm_add_application(cpchop SOURCES examples/cpchop.cpp)
opm_add_application(cpchop_depthtrend SOURCES examples/cpchop_depthtrend.cpp)
opm_add_application(cpregularize SOURCES examples/cpregularize.cpp)
opm_add_application(exp_variogram SOURCES examples/exp_variogram.cpp)
opm_add_application(grdecldips SOURCES examples/grdecldips.cpp)
opm_add_application(implicitcap_test SOURCES examples/implicitcap_test.cpp)
opm_add_application(known_answer_test SOURCES examples/known_answer_test.cpp)
opm_add_application(mimetic_aniso_solver_test SOURCES examples/mimetic_aniso_solver_test.cpp)
opm_add_application(mimetic_periodic_test SOURCES examples/mimetic_periodic_test.cpp)
opm_add_application(mimetic_solver_test SOURCES examples/mimetic_solver_test.cpp)
opm_add_application(sim_blackoil_impes SOURCES examples/sim_blackoil_impes.cpp)
opm_add_application(sim_co2_impes SOURCES examples/sim_co2_impes.cpp)
opm_add_application(sim_steadystate_explicit SOURCES examples/sim_steadystate_explicit.cpp)
opm_add_application(sim_steadystate_implicit SOURCES examples/sim_steadystate_implicit.cpp)
opm_add_application(steadystate_test_implicit SOURCES examples/steadystate_test_implicit.cpp)
opm_add_application(upscale_avg SOURCES examples/upscale_avg.cpp)
opm_add_application(upscale_cap SOURCES examples/upscale_cap.cpp)
opm_add_application(upscale_cond SOURCES examples/upscale_cond.cpp)
opm_add_application(upscale_perm SOURCES examples/upscale_perm.cpp)
opm_add_application(upscale_relperm SOURCES examples/upscale_relperm.cpp)
opm_add_application(upscale_relpermvisc SOURCES examples/upscale_relpermvisc.cpp)
opm_add_application(upscale_singlephase SOURCES examples/upscale_singlephase.cpp)
opm_add_application(upscale_steadystate_implicit SOURCES examples/upscale_steadystate_implicit.cpp)
opm_add_application(compareUpscaling SOURCES tests/compareUpscaling.cpp)
opm_add_application(upscale_relperm_benchmark SOURCES benchmarks/upscale_relperm_benchmark.cpp
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")
opm_add_application(upscale_elasticity SOURCES examples/upscale_elasticity.cpp
  LIBRARIES "${Boost_LIBRARIES}" INCLUDE_DIRS "${Boost_INCLUDE_DIRS}")

# setup extra tests (using helper binaries)
include (${CMAKE_CURRENT_SOURCE_DIR}/compareUpscaling.cmake)

# encode test cases so they can be embedded in the benchmark executables
include(EmbedCases.cmake)

# Setup static benchmarks
include(OpmStaticTargets)
opm_from_git(${PROJECT_SOURCE_DIR} benchmarks ${VCS_SHA1} benchmarks)
add_dependencies(benchmarks-static opm-grid-static)

# Copy static benchmarks to main project binary directory
foreach(benchmark ${OPM_BENCHMARKS})
  add_custom_command(TARGET benchmarks-static POST_BUILD
                     COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_BINARY_DIR}/static/benchmarks/src/benchmarks-static-build/bin/${benchmark}
                                                      ${CMAKE_BINARY_DIR}/bin/${benchmark}-static)
	if(INSTALL_BENCHMARKS)
		install(TARGETS ${benchmark} DESTINATION bin)
	endif()
endforeach()

# finalize the dune project, e.g. generating config.h etc.
finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
