# -*- mode: cmake; tab-width: 2; indent-tabs-mode: nil; truncate-lines: t; compile-command: "cmake -Wdev" -*-
# vim: set filetype=cmake autoindent tabstop=2 shiftwidth=2 expandtab softtabstop=2 nowrap:

###########################################################################
#                                                                         #
# Note: The bulk of the build system is located in the cmake/ directory.  #
#       This file only contains the specializations for this particular   #
#       project. Most likely you are interested in editing one of these   #
#       files instead:                                                    #
#                                                                         #
#       dune.module                              Name and version number  #
#       CMakeLists_files.cmake                   Path of source files     #
#       cmake/Modules/${project}-prereqs.cmake   Dependencies             #
#                                                                         #
###########################################################################

cmake_minimum_required (VERSION 3.23)

# Mandatory call to project
project(opm-upscaling C CXX)

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
option(INSTALL_BENCHMARKS "Install benchmark applications?" OFF)

macro(opm-upscaling_fortran_hook)
  #  interface requires a variable to inspect, a value won't suffice
  set(UNCONDITIONALLY TRUE)
  set(opm-upscaling_FORTRAN_IF UNCONDITIONALLY)
endmacro()

macro(opm-upscaling_config_hook)
  opm_need_version_of ("dune-common")
  opm_need_version_of ("dune-geometry")
  opm_need_version_of ("dune-istl")
endmacro()

macro(opm-upscaling_tests_hook)
  # setup extra tests (using helper binaries)
  include(compareUpscaling.cmake)
endmacro()

macro(opm-upscaling_targets_hook)
  # pack these cases which are alternatives in the code
  if(TARGET Boost::iostreams)
    include(EmbedCases.cmake)
    opm_pack_stone(upscale_relperm_benchmark)
    opm_pack_case(upscale_relperm_benchmark benchmark20)

    if(INSTALL_BENCHMARKS)
      add_custom_target(benchmarks ALL DEPENDS upscale_relperm_benchmark)
      set_target_properties(upscale_relperm_benchmark PROPERTIES EXCLUDE_FROM_ALL 0)
    else()
      add_custom_target(benchmarks DEPENDS upscale_relperm_benchmark)
    endif()
  endif()
endmacro()

macro(opm-upscaling_install_hook)
  install(
    DIRECTORY
      doc/man1
    DESTINATION
      ${CMAKE_INSTALL_MANDIR}
    FILES_MATCHING PATTERN
      "*.1"
  )
endmacro()

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
      set(opm-common_DIR "${_parent_full_dir}/opm-common")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

include(OpmInit)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)
OpmSetPolicies()

# list of prerequisites for this particular project; this is in a
# separate file (in cmake/Modules sub-directory) because it is shared
# with the find module
include (${project}-prereqs)

# read the list of components from this file (in the project directory);
# it should set various lists with the names of the files to include
include (CMakeLists_files.cmake)

# all setup common to the OPM library modules is done here
include (OpmLibMain)
