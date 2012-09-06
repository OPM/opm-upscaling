# -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
AC_PREREQ(2.50)
DUNE_AC_INIT # gets module version from dune.module file
AM_INIT_AUTOMAKE

# Use "silent" rules by default if available in the Automake being used to
# process this setup.  See, e.g.,
#
#   http://www.flameeyes.eu/autotools-mythbuster/automake/silent.html
#
# for details, particularly concerning the backwards compatibility.
#
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CONFIG_SRCDIR([dune/porsol/mimetic/IncompFlowSolverHybrid.hpp])
AM_CONFIG_HEADER([config.h])


# we need no more than the standard DE-stuff
# this module depends on dune-common dune-grid dune-istl dune-cornerpoint
# this implies checking for [dune-common], [dune-grid], [dune-istl], [dune-cornerpoint]
DUNE_CHECK_ALL

# Whether or not to enable features that depend on "constexpr" support.
AM_CONDITIONAL([HAVE_CXX0X_CONSTEXPR],dnl
  [test x"$HAVE_CXX0X_CONSTEXPR" = x"yes"])[]dnl

# implicitly set the Dune-flags everywhere
AC_SUBST(AM_CPPFLAGS, '$(DUNE_CPPFLAGS)  -I$(top_srcdir)')
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
  dune/porsol/blackoil/Makefile
  dune/porsol/blackoil/fluid/Makefile
  dune/porsol/blackoil/test/Makefile
  dune/porsol/twophase2/Makefile
  dune/porsol/twophase2/OPMIRISCode/Makefile
  dune/porsol/twophase2/OPMKvasiSophusCode/Makefile
  lib/Makefile
  examples/Makefile
  opm-porsol.pc
])

#  dune/porsol/opmtransport/Makefile
#  dune/porsol/opmtransport/src/Makefile
#  dune/porsol/opmtransport/examples/Makefile

AC_OUTPUT

# finally print the summary information
DUNE_SUMMARY_ALL
