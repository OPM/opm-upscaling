# $Date$
# $Revision$

check_PROGRAMS = 
noinst_PROGRAMS = upscale_perm \
                  upscale_relperm \
                  upscale_avg  \
                  upscale_cond \
                  cpchop \
                  exp_variogram \
                  cpregularize

upscale_perm_SOURCES = upscale_perm.cpp

upscale_relperm_SOURCES = upscale_relperm.cpp

upscale_cond_SOURCES = upscale_cond.cpp

upscale_avg_SOURCES = upscale_avg.cpp

cpchop_SOURCES = cpchop.cpp

exp_variogram_SOURCES = exp_variogram.cpp

cp_regularize_SOURCES = cp_regularize.cpp

AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNEMPILIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)   \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)


TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
