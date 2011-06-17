# $Date$
# $Revision$

check_PROGRAMS =
noinst_PROGRAMS = bo_fluid_test bo_well_test bo_fluid_pressuredeps bo_fluid_p_and_z_deps

bo_fluid_test_SOURCES = bo_fluid_test.cpp
bo_fluid_pressuredeps_SOURCES = bo_fluid_pressuredeps.cpp
bo_well_test_SOURCES = bo_well_test.cpp
bo_fluid_p_and_z_deps_SOURCES = bo_fluid_p_and_z_deps.cpp


AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(ALUGRID_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
        $(BOOST_SYSTEM_LIB) $(DUNEMPILIBS)              \
	 $(ALUGRID_LIBS) $(SUPERLU_LIBS)                \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)            \
	../fluid/libblackoilfluid_noinst.la


TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
