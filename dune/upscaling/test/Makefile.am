# $Date$
# $Revision$
#        aniso_implicit_steadystate_test \
#        aniso_steadystate_test \
#        implicit_steadystate_test \
#        opmpressure_test \
#        steadystate_test \
#        steadystate_test_implicit \
#        upscaling_test 

check_PROGRAMS = 
noinst_PROGRAMS = \
		 test_array \
         steadystate_test \ 
         steadystate_test_implicit
        

AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNEMPILIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)   \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)

test_array_SOURCES = test_array.cpp

#upscaling_test_SOURCES = upscaling_test.cpp

#opmpressure_test_SOURCES = opmpressure_test.cpp

steadystate_test_SOURCES = steadystate_test.cpp

steadystate_test_implicit_SOURCES = steadystate_test_implicit.cpp

#aniso_steadystate_test_SOURCES = aniso_steadystate_test.cpp

#implicit_steadystate_test_SOURCES = implicit_steadystate_test.cpp

#aniso_implicit_steadystate_test_SOURCES = aniso_implicit_steadystate_test.cpp


TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
