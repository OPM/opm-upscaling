# $Date$
# $Revision$

check_PROGRAMS =
noinst_PROGRAMS =				\
aniso_implicit_steadystate_test			\
aniso_steadystate_test				\
implicit_steadystate_test			\
opmpressure_test				\
steadystate_test_explicit			\
steadystate_test_implicit			\
upscaling_test

#test_array

AM_CPPFLAGS += $(DUNEMPICPPFLAGS) $(BOOST_CPPFLAGS) $(SUPERLU_CPPFLAGS)
AM_LDFLAGS  += $(DUNEMPILDFLAGS) $(BOOST_LDFLAGS) $(SUPERLU_LDFLAGS)
LDADD = $(DUNEMPILIBS) $(BOOST_UNIT_TEST_FRAMEWORK_LIB) \
        $(BOOST_FILESYSTEM_LIB) $(BOOST_SYSTEM_LIB)   \
        $(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS) $(SUPERLU_LIBS)


aniso_implicit_steadystate_test_SOURCES =	\
aniso_implicit_steadystate_test.cpp

aniso_steadystate_test_SOURCES =		\
aniso_steadystate_test.cpp

implicit_steadystate_test_SOURCES =		\
implicit_steadystate_test.cpp

opmpressure_test_SOURCES =			\
opmpressure_test.cpp

steadystate_test_explicit_SOURCES =			\
steadystate_test_explicit.cpp

steadystate_test_implicit_SOURCES =		\
steadystate_test_implicit.cpp

upscaling_test_SOURCES =			\
upscaling_test.cpp

#test_array_SOURCES = test_array.cpp

TESTS = $(check_PROGRAMS)

include $(top_srcdir)/am/global-rules
