# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([dune/porsol/mimetic/IncompFlowSolverHybrid.hpp])
AM_CONFIG_HEADER([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid dune-istl dune-cornerpoint
# this implies checking for [dune-common], [dune-grid], [dune-istl], [dune-cornerpoint]
DUNE_CHECK_ALL

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, $DUNE_CPPFLAGS)
AC_SUBST(AM_LDFLAGS, $DUNE_LDFLAGS)
LIBS="$DUNE_LIBS"

AC_CONFIG_FILES([
  Makefile
  doc/Makefile
  doc/doxygen/Makefile
  doc/doxygen/Doxyfile
  m4/Makefile
  dune/Makefile
  dune/porsol/Makefile
  dune/porsol/common/Makefile
  dune/porsol/common/test/Makefile
  dune/porsol/euler/Makefile
  dune/porsol/euler/test/Makefile
  dune/porsol/mimetic/Makefile
  dune/porsol/mimetic/test/Makefile
  examples/Makefile
  dune-porsol.pc
])
AC_OUTPUT
# finally print the summary information
DUNE_SUMMARY_ALL
