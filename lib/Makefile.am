# $Date$
# $Revision$

# The dune-porsol library
lib_LTLIBRARIES = libduneporsol.la

libduneporsol_la_SOURCES =

# Dummy C++ source to cause C++ linking.
nodist_EXTRA_libduneporsol_la_SOURCES = dummy.cpp

# Collect constituent parts into module library.
libduneporsol_la_LIBADD = \
   ../dune/porsol/opmpressure/src/libopmpressure_noinst.la \
   ../dune/porsol/common/libduneporsolcommon.la \
   ../dune/porsol/blackoil/fluid/libblackoilfluid_noinst.la \
   ../dune/porsol/twophase2/OPMKvasiSophusCode/libOPMKvasiSophusCode_noinst.la \
   ../dune/porsol/twophase2/OPMIRISCode/libOPMIRISCode_noinst.la

include $(top_srcdir)/am/global-rules
