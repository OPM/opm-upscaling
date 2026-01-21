# defines that must be present in config.h for our headers
set (opm-upscaling_CONFIG_VAR
  HAVE_SUPERLU
  HAVE_OPENMP
  )

# dependencies
set (opm-upscaling_DEPS
  # various runtime library enhancements
  "Boost 1.44.0
     COMPONENTS date_time system unit_test_framework REQUIRED"
  "Boost 1.44.0 COMPONENTS iostreams"
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
