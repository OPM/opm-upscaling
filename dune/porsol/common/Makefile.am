# $Date$
# $Revision$

SUBDIRS = test

commondir = $(includedir)/dune/porsol/common
common_HEADERS = blas_lapack.hpp BoundaryConditions.hpp fortran.hpp \
                 GridInterfaceEuler.hpp linearInterpolation.hpp Matrix.hpp \
                 MatrixInverse.hpp NonuniformTableLinear.hpp PeriodicHelpers.hpp \
                 ReservoirPropertyCapillary.hpp setupBoundaryConditions.hpp \
                 setupGridAndProps.hpp SimulatorBase.hpp \
		 SimulatorUtilities.hpp SintefLegacyGridInterface.hpp

noinst_LTLIBRARIES = libduneporsolcommon.la

libduneporsolcommon_la_SOURCES = \
	LinearSolverISTL.cpp


libduneporsolcommon_la_CPPFLAGS = $(DUNE_CPPFLAGS) $(BOOST_CPPFLAGS)
libduneporsolcommon_la_LDFLAGS  = $(DUNE_LDFLAGS) $(BOOST_LDFLAGS)

include $(top_srcdir)/am/global-rules
