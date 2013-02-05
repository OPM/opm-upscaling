# $Date$
# $Revision$

check_PROGRAMS = 
noinst_PROGRAMS =				\
cpchop						\
cpchop_depthtrend				\
cpregularize					\
exp_variogram					\
grdecldips                   			\
upscale_avg					\
upscale_cap					\
upscale_cond					\
upscale_elasticity				\
upscale_perm					\
upscale_relperm_ssor                            \
upscale_relperm_bgs                             \
upscale_relperm_aniso                           \
upscale_relperm_aniso_ssor                      \
upscale_relperm					\
upscale_relpermvisc                             \
upscale_steadystate_implicit


cpchop_SOURCES = cpchop.cpp

cpchop_depthtrend_SOURCES = cpchop_depthtrend.cpp

cpregularize_SOURCES = cpregularize.cpp

exp_variogram_SOURCES = exp_variogram.cpp

grdecldips_SOURCES = grdecldips.cpp

upscale_avg_SOURCES = upscale_avg.cpp

upscale_cap_SOURCES = upscale_cap.cpp

upscale_cond_SOURCES = upscale_cond.cpp

upscale_elasticity_SOURCES = upscale_elasticity.cpp
upscale_elasticity_CXXFLAGS = ${OPENMP_CXXFLAGS}
upscale_elasticity_LDADD = \
$(top_srcdir)/opm/elasticity/libelasticityupscale.la \
$(LDADD)

upscale_perm_SOURCES = upscale_perm.cpp

upscale_relperm_SOURCES = upscale_relperm.cpp

upscale_relperm_ssor_SOURCES = upscale_relperm.cpp
upscale_relperm_ssor_CXXFLAGS = -DSMOOTHER_ILU=0 $(AM_CPPFLAGS)

upscale_relperm_bgs_SOURCES = upscale_relperm.cpp
upscale_relperm_bgs_CXXFLAGS = -DSMOOTHER_ILU=0 -DSMOOTHER_BGS=1 $(AM_CPPFLAGS)

upscale_relperm_aniso_SOURCES = upscale_relperm.cpp
upscale_relperm_aniso_CXXFLAGS = -DANISOTROPIC_3D=1

upscale_relperm_aniso_ssor_SOURCES = upscale_relperm.cpp
upscale_relperm_aniso_ssor_CXXFLAGS = -DANISOTROPIC_3D=1 -DSMOOTHER_ILU=0 $(AM_CPPFLAGS)

upscale_relpermvisc_SOURCES = upscale_relpermvisc.cpp

upscale_steadystate_implicit_SOURCES = upscale_steadystate_implicit.cpp

AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNEMPILIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)   \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)


TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
