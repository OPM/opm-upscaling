# -*- mode: cmake; tab-width: 2; indent-tabs-mode: t; truncate-lines: t
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 noexpandtab softtabstop=2 nowrap:

cmake_minimum_required(VERSION 3.0)


# specify the option that allows to select the build system. if
# dunecontrol is used, we default to the dune build system. if not, we
# use the OPM home-brewn one
set(USE_DUNE_BUILDSYSTEM_DEFAULT OFF)
if (DEFINED ENV{DUNE_CONTROL_PATH})
  set(USE_DUNE_BUILDSYSTEM_DEFAULT ON)
endif()
option(USE_DUNE_BUILDSYSTEM "Use the Dune build system?" ${USE_DUNE_BUILDSYSTEM_DEFAULT})

if (USE_DUNE_BUILDSYSTEM)

######
# use DUNE's build system
######

# set up project and specify the minimum cmake version
project("opm-upscaling" C CXX)

option(INSTALL_BENCHMARKS "Install benchmark applications?" OFF)

# Sibling build
option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
if(SIBLING_SEARCH)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  foreach(MOD dune-common dune-geometry dune-grid dune-localfunctions dune-istl opm-common libecl opm-parser opm-grid opm-material opm-core)
    # do not overwrite the location of a module if it was explicitly
    # specified by the user
    if(${MOD}_DIR)
      continue()
    endif()

    # Try various possible locations for the build directory of the dependency
    foreach(BUILD_DIR "${_leaf_dir_name}" "build-cmake" "build" "../${MOD}-build" "../build-${MOD}" "../build" ".")
      if(EXISTS "${_modules_dir}/${MOD}/${BUILD_DIR}/${MOD}-config.cmake")
        set(${MOD}_DIR ${_modules_dir}/${MOD}/${BUILD_DIR})
        break()
      endif()
    endforeach()
  endforeach()
endif()

if(dune-common_DIR AND NOT IS_DIRECTORY ${dune-common_DIR})
  message(WARNING "Value ${dune-common_DIR} passed to variable"
    " dune-common_DIR is not a directory")
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

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

else()

######
# use the traditional OPM build system
######

# Mandatory call to project
project(opm-upscaling CXX)

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
option(INSTALL_BENCHMARKS "Install benchmark applications?" OFF)

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common}")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

include(OpmInit)

# not the same location as most of the other projects? this hook overrides
macro (dir_hook)
endmacro (dir_hook)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

set(USE_OPENMP_DEFAULT ON)

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include (${project}-prereqs)

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

macro (config_hook)
	opm_need_version_of ("dune-common")
	opm_need_version_of ("dune-geometry")
	opm_need_version_of ("dune-istl")
	list (APPEND ${project}_CONFIG_IMPL_VARS
		HAVE_LAPACK
		HAVE_UMFPACK
		HAVE_SUPERLU
		HAVE_SUITESPARSE_UMFPACK
		)
endmacro (config_hook)

macro (prereqs_hook)
endmacro (prereqs_hook)

macro (sources_hook)
endmacro (sources_hook)

macro (files_hook)
endmacro (files_hook)

macro (fortran_hook)
	#  interface requires a variable to inspect, a value won't suffice
	set (UNCONDITIONALLY TRUE)
	set (${project}_FORTRAN_IF UNCONDITIONALLY)
endmacro (fortran_hook)

macro (tests_hook)
endmacro (tests_hook)

# all setup common to the OPM library modules is done here
include (OpmLibMain)

# setup extra tests (using helper binaries)
include (${CMAKE_CURRENT_SOURCE_DIR}/compareUpscaling.cmake)

# encode test cases so they can be embedded in the benchmark executables
include (${PROJECT_SOURCE_DIR}/EmbedCases.cmake)

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

endif()
