SUBDIRS =

check_PROGRAMS =
noinst_PROGRAMS =				\
aniso_implicitcap_test				\
aniso_simulator_test				\
blackoil_sim_test				\
grdecl_to_legacy_test				\
ifsh_solver_test				\
implicitcap_test				\
istl_test					\
known_answer_test				\
make_vtk_test					\
mimetic_aniso_solver_test			\
mimetic_periodic_test				\
mimetic_solver_test				\
simulator_test					\
simulator_implicit_test				\
spe10_test					\
spu_2p						\
tpfa_compressible_solver_test			\
tpfa_solver_test
#				\
#twophase2_test



grdecl_to_legacy_test_SOURCES = grdecl_to_legacy_test.cpp

make_vtk_test_SOURCES  = make_vtk_test.cpp

ifsh_solver_test_SOURCES = ifsh_solver_test.cpp
ifsh_solver_test_LDADD = $(LDADD) ../lib/libopmporsol.la

mimetic_solver_test_SOURCES = mimetic_solver_test.cpp

mimetic_aniso_solver_test_SOURCES = mimetic_aniso_solver_test.cpp

mimetic_periodic_test_SOURCES = mimetic_periodic_test.cpp

istl_test_SOURCES = istl_test.cpp

spu_2p_SOURCES = spu_2p.cpp
spu_2p_LDADD = $(LDADD) ../lib/libopmporsol.la

spe10_test_SOURCES = spe10_test.cpp

tpfa_solver_test_SOURCES = tpfa_solver_test.cpp
tpfa_solver_test_LDADD = $(LDADD) ../lib/libopmporsol.la

tpfa_compressible_solver_test_SOURCES = tpfa_compressible_solver_test.cpp
tpfa_compressible_solver_test_LDADD = $(LDADD) ../lib/libopmporsol.la

known_answer_test_SOURCES = known_answer_test.cpp

simulator_test_SOURCES = simulator_test.cpp
#simulator_test_LDADD   = $(LDADD) -ltbb

simulator_implicit_test_SOURCES = simulator_implicit_test.cpp

aniso_simulator_test_SOURCES = aniso_simulator_test.cpp

implicitcap_test_SOURCES = implicitcap_test.cpp

aniso_implicitcap_test_SOURCES = aniso_implicitcap_test.cpp

blackoil_sim_test_SOURCES = blackoil_sim_test.cpp
blackoil_sim_test_LDADD = $(LDADD) ../lib/libopmporsol.la

#twophase2_test_SOURCES = twophase2_test.cpp
#twophase2_test_LDADD = $(LDADD) ../lib/libopmporsol.la

# ----------------------------------------------------------------------
# Tests programs reliant on recent compiler support, notably the
# "constexpr" feature.
#
if HAVE_CXX0X_CONSTEXPR
noinst_PROGRAMS +=				\
co2_sim_test                                    \
co2_blackoil_pvt

co2_sim_test_SOURCES = co2_sim_test.cpp
co2_sim_test_LDADD = $(LDADD) ../lib/libopmporsol.la
co2_blackoil_pvt_SOURCES = co2_blackoil_pvt.cpp
co2_blackoil_pvt_LDADD = $(LDADD) ../lib/libopmporsol.la
endif


TESTS = $(check_PROGRAMS)

# Commented out ALUGRID flags and libs, to avoid link trouble
# related to debug vs. non-debug versions of std::vector used
# in the ALUgrid interface.
# If we find out how to properly install two versions of
# ALUgrid (both debug and release), we may reinstate them.
# AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(OPM_BOOST_CPPFLAGS) $(ALUGRID_CPPFLAGS) $(SUPERLU_CPPFLAGS)
# AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(OPM_BOOST_LDFLAGS) $(ALUGRID_LDFLAGS) $(SUPERLU_LDFLAGS)
# LDADD = $(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
          $(BOOST_SYSTEM_LIB) $(DUNEMPILIBS)            \
          $(ALUGRID_LIBS) $(SUPERLU_LIBS)               \
          $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(OPM_BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(OPM_BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)

LDADD =                                         \
$(DUNE_LIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB)   \
$(BOOST_SYSTEM_LIB) $(DUNEMPILIBS)              \
$(SUPERLU_LIBS)                                 \
$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)

# don't follow the full GNU-standard
# we need automake 1.5
AUTOMAKE_OPTIONS = foreign 1.5

# pass most important options when "make distcheck" is used
DISTCHECK_CONFIGURE_FLAGS =                             \
--with-dune-common=$(DUNE_COMMON_ROOT)                  \
--with-dune-grid=$(DUNE_GRID_ROOT)                      \
--with-dune-istl=$(DUNE_ISTL_ROOT)                      \
--with-dune-cornerpoint=$(DUNE_CORNERPOINT_ROOT)        \
CXX="$(CXX)" CC="$(CC)"

include $(top_srcdir)/am/global-rules
