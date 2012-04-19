AC_DEFUN([AX_CXX0X_CONSTEXPR],
[
AC_REQUIRE([GXX0X])

AC_CACHE_CHECK([for constexpr],
  [ax_cv_have_constexpr],
  [AS_IF([test x"$HAVE_CXX0X" != x"no"],
    [
     AC_LANG_PUSH([C++])[]dnl

     AC_LINK_IFELSE(
       AC_LANG_PROGRAM([[
         struct S {
           static double constexpr s = 1.0;
         };]],dnl
         [[
           S s;
         ]]),dnl
         [ax_cv_have_constexpr="yes"],dnl
         [ax_cv_have_constexpr="no"]
       )[]dnl

     AC_LANG_POP([C++])[]dnl
    ],dnl
    [ax_cv_have_constexpr="no"]
    )[]dnl
  ])[]dnl

HAVE_CXX0X_CONSTEXPR="$ax_cv_have_constexpr"
])[]dnl
