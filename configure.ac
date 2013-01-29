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

AC_CONFIG_SRCDIR([opm/porsol/mimetic/IncompFlowSolverHybrid.hpp])
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
  opm/Makefile
  opm/porsol/Makefile
  opm/porsol/common/Makefile
  opm/porsol/common/test/Makefile
  opm/porsol/euler/Makefile
  opm/porsol/euler/test/Makefile
  opm/porsol/mimetic/Makefile
  opm/porsol/mimetic/test/Makefile
  opm/porsol/blackoil/Makefile
  opm/porsol/blackoil/fluid/Makefile
  opm/porsol/blackoil/test/Makefile
  opm/porsol/twophase2/Makefile
  opm/porsol/twophase2/OPMIRISCode/Makefile
  opm/porsol/twophase2/OPMKvasiSophusCode/Makefile
  lib/Makefile
  examples/Makefile
  opm-porsol.pc
])

#  opm/porsol/opmtransport/Makefile
#  opm/porsol/opmtransport/src/Makefile
#  opm/porsol/opmtransport/examples/Makefile

AC_OUTPUT

# finally print the summary information
DUNE_SUMMARY_ALL
