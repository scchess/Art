AC_INIT([art C++], [GSM.04.18.2016], [whduke@gmail.com], [art], [http://http://www.niehs.nih.gov/research/resources/software/biostatistics/art])
AC_PREREQ([2.69])
#AC_PREREQ([2.63])
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AM_INIT_AUTOMAKE([subdir-objects])
AC_PROG_CXX
AC_PROG_CC
# Checks for libraries.
# Checks for header files.
AC_CHECK_HEADERS([stdlib.h])
AC_CHECK_LIB(m,main)
AC_CHECK_LIB(gslcblas,main,,AC_MSG_ERROR([gslcblas not found]))
AC_CHECK_LIB(gsl,main,,AC_MSG_ERROR([gsl library not found]))
AC_PROG_MAKE_SET

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

