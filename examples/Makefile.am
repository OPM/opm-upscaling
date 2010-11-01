SUBDIRS =

check_PROGRAMS =
noinst_PROGRAMS = grdecl_to_legacy_test \
                  make_vtk_test \
                  istl_test \
                  ifsh_solver_test \
                  mimetic_solver_test \
                  mimetic_aniso_solver_test \
                  mimetic_periodic_test \
		  spe10_test \
		  tpfa_solver_test \
		  known_answer_test

grdecl_to_legacy_test_SOURCES = grdecl_to_legacy_test.cpp

make_vtk_test_SOURCES  = make_vtk_test.cpp

ifsh_solver_test_SOURCES = ifsh_solver_test.cpp
ifsh_solver_test_LDADD = $(LDADD) ../lib/libduneporsol.la

mimetic_solver_test_SOURCES = mimetic_solver_test.cpp

mimetic_aniso_solver_test_SOURCES = mimetic_aniso_solver_test.cpp

mimetic_periodic_test_SOURCES = mimetic_periodic_test.cpp

istl_test_SOURCES = istl_test.cpp

spe10_test_SOURCES = spe10_test.cpp

tpfa_solver_test_SOURCES = tpfa_solver_test.cpp
tpfa_solver_test_LDADD = $(LDADD) ../lib/libduneporsol.la

known_answer_test_SOURCES = known_answer_test.cpp

TESTS = $(check_PROGRAMS)


AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(ALUGRID_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
        $(BOOST_SYSTEM_LIB) $(DUNEMPILIBS)              \
	 $(ALUGRID_LIBS) $(SUPERLU_LIBS)                \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS = --with-dune-common=$(DUNE_COMMON_ROOT) \
                            --with-dune-grid=$(DUNE_GRID_ROOT) \
                            --with-dune-istl=$(DUNE_ISTL_ROOT) \
                            --with-dune-cornerpoint=$(DUNE_CORNERPOINT_ROOT) \
                            CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
