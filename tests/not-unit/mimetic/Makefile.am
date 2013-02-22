# $Date$
# $Revision$

check_PROGRAMS =
noinst_PROGRAMS = mimetic_ipeval_test

AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(OPM_BOOST_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(OPM_BOOST_LDFLAGS) $(ALUGRID_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
        $(BOOST_SYSTEM_LIB) $(DUNEMPILIBS)              \
	 $(ALUGRID_LIBS) $(SUPERLU_LIBS)                \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

mimetic_ipeval_test_SOURCES = mimetic_ipeval_test.cpp


TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
