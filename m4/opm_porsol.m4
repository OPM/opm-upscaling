dnl -*- autoconf -*-

dnl locate opm-porsol library itself; this macro is called by every module
dnl that depends on opm-porsol.
AC_DEFUN([OPM_PORSOL_CHECK_MODULE],
[
 OPM_CHECK_PKG_MODULE([opm-porsol],[1.0],[DUNE module containing porous media PDE solvers])
])

dnl find all prerequisites of opm-porsol; nothing to do here since this
dnl is done by the CMake module and then stored in the -config file.
AC_DEFUN([OPM_PORSOL_CHECKS],[])
