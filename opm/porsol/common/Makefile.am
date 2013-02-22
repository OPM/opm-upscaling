# $Date$
# $Revision$

SUBDIRS =

commondir = $(includedir)/dune/porsol/common
common_HEADERS = blas_lapack.hpp BoundaryConditions.hpp fortran.hpp \
                 GridInterfaceEuler.hpp Matrix.hpp \
                 MatrixInverse.hpp NonuniformTableLinear.hpp PeriodicHelpers.hpp \
                 ReservoirPropertyCapillary.hpp setupBoundaryConditions.hpp \
                 setupGridAndProps.hpp SimulatorBase.hpp \
		 SimulatorUtilities.hpp SintefLegacyGridInterface.hpp

noinst_LTLIBRARIES = libopmporsolcommon.la

libopmporsolcommon_la_SOURCES = \
	LinearSolverISTL.cpp


libopmporsolcommon_la_CPPFLAGS = $(AM_CPPFLAGS) $(DUNE_CPPFLAGS) $(OPM_BOOST_CPPFLAGS)
libopmporsolcommon_la_LDFLAGS  = $(DUNE_LDFLAGS) $(OPM_BOOST_LDFLAGS)

include $(top_srcdir)/am/global-rules
