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
# opm-upscaling has no Fortran sources; Fortran is enabled only so that
# FortranCInterface can detect the BLAS/LAPACK name mangling and write
# FCMacros.h, which opm/porsol/common/blas_lapack.hpp hard-requires.
#
# Ask for the language only when a compiler for it actually exists, rather than
# keying on the host OS: a Windows machine with Intel Fortran should get the
# real header, and a Linux machine without gfortran should still configure.
# opm-upscaling_language_hook writes a fallback FCMacros.h when it is absent.
# Note: project() has not run yet, so use check_language(), not CMAKE_Fortran_COMPILER.
include(CheckLanguage)
check_language(Fortran)
if(CMAKE_Fortran_COMPILER)
  project(opm-upscaling C CXX Fortran)
else()
  project(opm-upscaling C CXX)
endif()

option(SIBLING_SEARCH "Search for other modules in sibling directories?" ON)
option(INSTALL_BENCHMARKS "Install benchmark applications?" OFF)

macro(opm-upscaling_language_hook)
  if(CMAKE_Fortran_COMPILER)
    include(FortranCInterface)
    fortrancinterface_header(${PROJECT_BINARY_DIR}/FCMacros.h MACRO_NAMESPACE FC_)
  else()
    # Without a Fortran compiler FortranCInterface cannot detect the BLAS/LAPACK
    # name mangling, but opm/porsol/common/blas_lapack.hpp hard-requires
    # FCMacros.h (#error when it is missing). Generate a fallback that assumes
    # the de-facto standard convention -- lowercase with a trailing underscore --
    # used by gfortran-built BLAS/LAPACK (OpenBLAS, reference LAPACK) and also
    # exported by MKL.
    # Only write the fallback if the file does not exist yet, so that a
    # hand-edited FCMacros.h (for a library with a different convention)
    # survives reconfiguration.
    if(NOT EXISTS ${PROJECT_BINARY_DIR}/FCMacros.h)
      message(STATUS "No Fortran compiler available - writing fallback FCMacros.h "
                     "assuming lowercase-with-underscore BLAS/LAPACK name mangling "
                     "(gfortran/OpenBLAS/LAPACK/MKL convention). If linking BLAS/LAPACK "
                     "fails with unresolved symbols like 'dgemm_', your library uses a "
                     "different convention; edit ${PROJECT_BINARY_DIR}/FCMacros.h to match.")
      file(WRITE ${PROJECT_BINARY_DIR}/FCMacros.h
"#ifndef FC_HEADER_INCLUDED
#define FC_HEADER_INCLUDED

/* Fallback Fortran name mangling, generated without a Fortran compiler:
   lowercase with trailing underscore (gfortran/OpenBLAS/LAPACK default). */
#define FC_GLOBAL(name,NAME) name##_
#define FC_GLOBAL_(name,NAME) name##_

#endif
")
    endif()
  endif()
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

if(SIBLING_SEARCH AND NOT opm-common_DIR)
  # guess the sibling dir
  get_filename_component(_leaf_dir_name ${PROJECT_BINARY_DIR} NAME)
  get_filename_component(_parent_full_dir ${PROJECT_BINARY_DIR} DIRECTORY)
  get_filename_component(_parent_dir_name ${_parent_full_dir} NAME)
  #Try if <module-name>/<build-dir> is used
  get_filename_component(_modules_dir ${_parent_full_dir} DIRECTORY)
  if(IS_DIRECTORY ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-common_DIR ${_modules_dir}/opm-common/${_leaf_dir_name})
    set(opm-grid_DIR ${_modules_dir}/opm-grid/${_leaf_dir_name})
  else()
    string(REPLACE ${PROJECT_NAME} opm-common _opm_common_leaf ${_leaf_dir_name})
    if(NOT _leaf_dir_name STREQUAL _opm_common_leaf
        AND IS_DIRECTORY ${_parent_full_dir}/${_opm_common_leaf})
      # We are using build directories named <prefix><module-name><postfix>
      set(opm-common_DIR ${_parent_full_dir}/${_opm_common_leaf})
      string(REPLACE ${PROJECT_NAME} opm-grid _opm_grid_leaf ${_leaf_dir_name})
      set(opm-grid_DIR ${_parent_full_dir}/${_opm_grid_leaf})
    elseif(IS_DIRECTORY ${_parent_full_dir}/opm-common)
      # All modules are in a common build dir
      set(opm-common_DIR "${_parent_full_dir}/opm-common")
      set(opm-grid_DIR "${_parent_full_dir}/opm-grid")
    endif()
  endif()
endif()
if(opm-common_DIR AND NOT IS_DIRECTORY ${opm-common_DIR})
  message(WARNING "Value ${opm-common_DIR} passed to variable"
    " opm-common_DIR is not a directory")
endif()

find_package(opm-common REQUIRED)

# project information is in dune.module. Read this file and set variables.
# we cannot generate dune.module since it is read by dunecontrol before
# the build starts, so it makes sense to keep the data there then.
include (OpmInit)

# all setup common to the OPM library modules is done here
include (OpmLibMain)
