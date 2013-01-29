# $Date$
# $Revision$

check_PROGRAMS = matrix_test boundaryconditions_test nonuniformtablelinear_test uniformtablelinear_test
noinst_PROGRAMS = \
        gie_test \
        periodic_test \
        rockjfunc_test

# Common flags
AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(OPM_BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(OPM_BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
	$(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)     \
	$(SUPERLU_LIBS) $(DUNEMPILIBS)

# Additional quirks for individual programs
matrix_test_SOURCES = matrix_test.cpp

boundaryconditions_test_SOURCES = boundaryconditions_test.cpp

nonuniformtablelinear_test_SOURCES = nonuniformtablelinear_test.cpp

uniformtablelinear_test_SOURCES = uniformtablelinear_test.cpp

periodic_test_SOURCES = periodic_test.cpp

gie_test_SOURCES = gie_test.cpp

rockjfunc_test_SOURCES = rockjfunc_test.cpp
rockjfunc_test_LDADD   = $(LDADD) $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
