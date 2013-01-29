# $Date$
# $Revision$

# The opm-porsol library
lib_LTLIBRARIES = libopmporsol.la

libopmporsol_la_SOURCES =

# Dummy C++ source to cause C++ linking.
nodist_EXTRA_libopmporsol_la_SOURCES = dummy.cpp

# Collect constituent parts into module library.
libopmporsol_la_LIBADD =                                                   \
../opm/porsol/common/libopmporsolcommon.la                                \
../opm/porsol/blackoil/fluid/libblackoilfluid_noinst.la
#                    \
#../opm/porsol/twophase2/OPMKvasiSophusCode/libOPMKvasiSophusCode_noinst.la \
#../opm/porsol/twophase2/OPMIRISCode/libOPMIRISCode_noinst.la

include $(top_srcdir)/am/global-rules
