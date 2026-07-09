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
project(opm-upscaling C CXX Fortran)

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
option(INSTALL_BENCHMARKS "Install benchmark applications?" OFF)

macro(opm-upscaling_language_hook)
  include(FortranCInterface)
  fortrancinterface_header(${PROJECT_BINARY_DIR}/FCMacros.h MACRO_NAMESPACE FC_)
endmacro()

macro(opm-upscaling_prereqs_hook)
  if(NOT TARGET SuperLU::SuperLU AND NOT TARGET SuiteSparse::UMFPACK)
    message(FATAL_ERROR "Enable either SuperLU or UMFPACK")
  endif()

  target_link_libraries(opmupscaling
    PUBLIC
      opmcommon
      opmgrid
      Boost::date_time
      BLAS::BLAS
      LAPACK::LAPACK
      dunecommon
      dunegeometry
      dunegrid
      Dune::ISTL
  )
  if(MPI_FOUND)
    target_link_libraries(opmupscaling PUBLIC MPI::MPI_C MPI::MPI_CXX)
  endif()
  if(TARGET SuperLU::SuperLU)
    target_link_libraries(opmupscaling PUBLIC SuperLU::SuperLU)
  endif()
  if(TARGET SuiteSparse::UMFPACK)
    target_link_libraries(opmupscaling PUBLIC SuiteSparse::UMFPACK)
  endif()
endmacro()

macro(opm-upscaling_tests_hook)
  # setup extra tests (using helper binaries)
  include(compareUpscaling.cmake)
endmacro()

macro(opm-upscaling_targets_hook)
  target_include_directories(opmupscaling PRIVATE ${PROJECT_BINARY_DIR})

  # pack these cases which are alternatives in the code
  # cannot build boost::iostreams 1.66 as c++-20.
  if(TARGET Boost::iostreams AND Boost_VERSION VERSION_GREATER 1.66)
    include(EmbedCases.cmake)
    opm_pack_stone(upscale_relperm_benchmark)
    opm_pack_case(upscale_relperm_benchmark benchmark20)
    target_link_libraries(upscale_relperm_benchmark PRIVATE Boost::iostreams)

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

# These modules are special because they search with the first
# find_package call in here
set(initial_depends opm-common dune-common)

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  foreach(module ${initial_depends})
    if(IS_DIRECTORY ${_modules_dir}/${module}/${_leaf_dir_name})
      set(${module}_DIR ${_modules_dir}/${module}/${_leaf_dir_name})
    else()
      string(REPLACE ${PROJECT_NAME} ${module} _module_leaf ${_leaf_dir_name})
      if(NOT _leaf_dir_name STREQUAL _module_leaf
          AND IS_DIRECTORY ${_parent_full_dir}/${_module_leaf})
        # We are using build directories named <prefix><module-name><postfix>
        set(${module}_DIR ${_parent_full_dir}/${_module_leaf})
      elseif(IS_DIRECTORY ${_parent_full_dir}/${module})
        # All modules are in a common build dir
        set(${module}_DIR "${_parent_full_dir}/${module}")
      endif()
    endif()
  endforeach()
endif()
foreach(module ${initial_depends})
  if(${module}_DIR AND NOT IS_DIRECTORY ${${module}_DIR})
    message(WARNING "Value ${${module}_DIR} passed to variable"
      " ${module}_DIR is not a directory")
  endif()
endforeach()

find_package(opm-common REQUIRED)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# all setup common to the OPM library modules is done here
include (OpmLibMain)
