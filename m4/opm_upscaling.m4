dnl -*- autoconf -*-

dnl locate opm-upscaling library itself; this macro is called by every module
dnl that depends on opm-upscaling.
AC_DEFUN([OPM_UPSCALING_CHECK_MODULE],
[
 OPM_CHECK_PKG_MODULE([opm-upscaling],[1.0],[DUNE module containing single-phase and steady-state upscaling methods])
])

dnl find all prerequisites of opm-upscaling; nothing to do here since this
dnl is done by the CMake module and then stored in the -config file.
AC_DEFUN([OPM_UPSCALING_CHECKS],[])
