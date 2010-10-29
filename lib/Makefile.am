# $Date$
# $Revision$

# The dune-porsol library
lib_LTLIBRARIES = libduneporsol.la

libduneporsol_la_SOURCES =

# Dummy C++ source to cause C++ linking.
nodist_EXTRA_libduneporsol_la_SOURCES = dummy.cpp

# Collect constituent parts into cpg support library.
libduneporsol_la_LIBADD = \
   ../dune/porsol/opmpressure/libopmpressure_noinst.la

include $(top_srcdir)/am/global-rules
