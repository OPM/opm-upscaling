# defines that must be present in config.h for our headers
set (opm-upscaling_CONFIG_VAR
  HAVE_SUPERLU
  HAVE_OPENMP
  HAVE_LAPACK
  HAVE_SUITESPARSE_UMFPACK
)

# CMake 3.30.0 requires to find Boost in CONFIG mode
if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.30.0)
  set(_Boost_CONFIG_MODE CONFIG)
endif()

# dependencies
set (opm-upscaling_DEPS
  # various runtime library enhancements
  "Boost 1.44.0
    COMPONENTS date_time iostreams REQUIRED ${_Boost_CONFIG_MODE}"
  # matrix library
  "BLAS REQUIRED"
  "LAPACK REQUIRED"
  # solver
  "SuperLU"
  # DUNE dependency
  "dune-common REQUIRED"
  "dune-istl REQUIRED"
  "dune-geometry REQUIRED"
  "dune-grid REQUIRED"
  "opm-common REQUIRED"
  "opm-grid REQUIRED"
  )

find_package_deps(opm-upscaling)

if(NOT HAVE_ECL_INPUT OR NOT HAVE_ECL_OUTPUT)
  message(FATAL_ERROR "Eclipse input/output support required in opm-common")
endif()
