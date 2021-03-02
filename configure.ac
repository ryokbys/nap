#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([nap], [2.0.0], [kobayashi.ryo@nitech.ac.jp])
AC_CONFIG_SRCDIR([pmd/parallel_md.F])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_FC(mpiifort mpifrtpx mpifrt openmpif90 mpif90)
AC_PROG_CPP
AC_PATH_PROG([PYTHON],[python])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([pmd/makefile mkconf/makefile fitpot/run_pmd.sh
                 fitpot/run_smd.sh qmcl/makefile neb/makefile
                 fitpot/makefile nappy/pmd/makefile])
AC_OUTPUT
