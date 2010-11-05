# $Date$
# $Revision$

check_PROGRAMS =
noinst_PROGRAMS = bo_transport_test

bo_transport_test_SOURCES = bo_transport_test.cpp


AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(ALUGRID_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
        $(BOOST_SYSTEM_LIB) $(DUNEMPILIBS)              \
	 $(ALUGRID_LIBS) $(SUPERLU_LIBS)                \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)


TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
