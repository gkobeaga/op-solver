AC_PREREQ([2.69])
AC_INIT([solver],[0.1],[gorka@kobeaga.eus])

AC_CONFIG_SRCDIR([src/op-solver.h])

AC_CONFIG_MACRO_DIR([m4])

dnl AM_INIT_AUTOMAKE
AM_INIT_AUTOMAKE([subdir-objects])
AM_PROG_AR

AC_CONFIG_HEADERS([config.h])

dnl Check for programs
dnl CFLAGS="-O3"
AC_PROG_CC
AC_PROG_INSTALL
LT_INIT

dnl Checks for libraries.
AC_CONFIG_MACRO_DIR([m4])

LT_INIT

dnl Checks for header files.
AC_CHECK_HEADERS([ stdlib.h math.h string.h strings.h errno.h assert.h stddef.h unistd.h malloc.h sys/types.h sys/times.h ])

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debug code [default=no]])],
                              [debug="$withval"], [debug=no])

dnl CPLEX path
AC_ARG_VAR([CPLEX_PATH], [path to CPLEX [default=/opt/]])
if test -z "$CPLEX_PATH"; then
  if test -d "/opt/ibm"; then
    CPLEX_PATH='/opt/ibm/ILOG/CPLEX_Studio125/cplex'
  else
    AC_MSG_ERROR([CPLEX path not found.])
    fi
fi

dnl GMP path
#AC_ARG_VAR([GMP_PATH], [path to GMP [default=/opt/]])
#if test -z "$GMP_PATH"; then
#  if test -d "/opt/gmp-6.1.2"; then
#    GMP_PATH='/opt/gmp-6.1.2'
#  else
#    AC_MSG_ERROR([GMP path not found.])
#    fi
#fi

dnl Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_PID_T
AC_TYPE_SIZE_T

dnl Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile
                 src/lib/Makefile
                 src/prob/Makefile
                 src/prob/lib/Makefile
                 src/prob/lib/graph/Makefile
                 src/prob/lib/data/Makefile
                 src/prob/lib/data/nearest/kdtree/Makefile
                 src/prob/lp/Makefile
                 src/prob/ip/Makefile
                 src/prob/ip/exact/Makefile
                 src/prob/ip/exact/bac/Makefile
                 src/prob/kp/Makefile
                 src/prob/cp/Makefile
                 src/prob/cp/init/Makefile
                 src/prob/cp/heur/Makefile
                 src/prob/cp/heur/ea/Makefile
                 src/prob/cp/exact/Makefile
                 src/prob/cp/exact/bac/Makefile
                 src/prob/tsp/Makefile
                 src/prob/tsp/init/Makefile
                 src/prob/tsp/heur/Makefile
                 src/prob/op/Makefile
                 src/prob/op/init/Makefile
                 src/prob/op/heur/Makefile
                 src/prob/op/heur/ea/Makefile
                 src/prob/op/exact/Makefile
                 src/prob/op/exact/bac/Makefile
                 ])
AC_OUTPUT

AC_MSG_NOTICE([summary of build options:
               Version:           ${PACKAGE_STRING}
               Debug Build:       ${debug}
               Host type:         ${host}
               Install prefix:    ${prefix}
               Compiler:          ${CC} ${CFLAGS}
               CPLEX path:        ${CPLEX_PATH}
               GMP path:          ${GMP_PATH}
               Static libraries:  ${enable_static}
               Shared libraries:  ${enable_shared}
               Linker:            ${LD} ${LDFLAGS} ${LIBS}])
