AC_PREREQ([2.69])
AC_INIT([op-solver],[0.1],[gorka@kobeaga.eus])

AC_CONFIG_SRCDIR([src/op-solver.h])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([subdir-objects])
AM_PROG_AR

AC_CONFIG_HEADERS([config.h])

dnl CFLAGS="-O3"
AC_PROG_CC
AC_PROG_CC_STDC
AC_PROG_INSTALL
LT_INIT

dnl Checks for header files.
AC_CHECK_HEADERS([ stdlib.h math.h string.h strings.h errno.h assert.h stddef.h unistd.h malloc.h sys/types.h sys/times.h gmp.h])

AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debug code [default=no]])],
                              [debug="$withval"], [debug=no])

AC_ARG_ENABLE([dummy-lp], AS_HELP_STRING([--enable-dummy-lp], [Use Dummy LP Solver]),[],[enable_dummy_lp=no] )
AS_IF([test "x$enable_dummy_lp" != "xno"],
      [
       have_dummy_lp=yes
       AM_CONDITIONAL(HAVE_LP_SOLVER, true)
       AM_CONDITIONAL(HAVE_LP_SOLVER_DUMMY, true)
       AM_CONDITIONAL(HAVE_LP_SOLVER_CPLEX, false)
       ],
       [
        have_dummy_lp=no
        AM_CONDITIONAL(HAVE_LP_SOLVER_DUMMY, false)
        ])

AS_IF([test "x$have_dummy_lp" = "xyes"],
      [
       AC_SUBST([LP_SOLVER_API])
       LP_SOLVER_API=lib/dummy.c
       AC_SUBST(LP_SOLVER_LIB, [])
       AC_SUBST(LP_SOLVER_INCLUDE, [])
       ],
       [AS_IF([test "x$enable_dummy_lp" = "xyes"],
              [AC_MSG_ERROR([Dummy LP requested but not found])
              ])
       ])

AC_ARG_WITH([cplex], AS_HELP_STRING([--with-cplex], [Use CPLEX]),[],[with_cplex=no] )
AS_IF([test "x$with_cplex" != "xno"],
      [
       have_cplex=yes

       AM_CONDITIONAL(HAVE_LP_SOLVER, true)
       AM_CONDITIONAL(HAVE_LP_SOLVER_CPLEX, true)
       AM_CONDITIONAL(HAVE_LP_SOLVER_DUMMY, false)
       ],
       [
        have_cplex=no
        AS_IF([test "x$have_dummy_lp" = "xno"],
              [
               AM_CONDITIONAL(HAVE_LP_SOLVER, false)
               ],
               [
                AM_CONDITIONAL(HAVE_LP_SOLVER, true)
                ])
        AM_CONDITIONAL(HAVE_LP_SOLVER_CPLEX, false)
        ])

AS_IF([test "x$have_cplex" = "xyes"],
      [
       AC_SUBST(LP_SOLVER_LIB, [$with_cplex/lib/x86-64_linux/static_pic/libcplex.a])
       AC_SUBST(LP_SOLVER_INCLUDE, [$with_cplex/include])
       ],
       [AS_IF([test "x$with_cplex" = "xyes"],
              [AC_MSG_ERROR([cplex requested but not found])
              ])
       ])

AC_TYPE_PID_T
AC_TYPE_SIZE_T

AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([memset pow sqrt])

AC_CONFIG_FILES([Makefile src/Makefile
                 src/lib/Makefile
                 src/graph/Makefile
                 src/data/Makefile
                 src/data/nearest/kdtree/Makefile
                 src/lp/Makefile
                 src/ip/Makefile
                 src/ip/exact/Makefile
                 src/ip/exact/bac/Makefile
                 src/kp/Makefile
                 src/cp/Makefile
                 src/cp/init/Makefile
                 src/cp/heur/Makefile
                 src/cp/heur/ea/Makefile
                 src/cp/exact/Makefile
                 src/cp/exact/bac/Makefile
                 src/tsp/Makefile
                 src/tsp/init/Makefile
                 src/tsp/heur/Makefile
                 src/tsp/heur/linkern/Makefile
                 src/op/Makefile
                 src/op/init/Makefile
                 src/op/heur/Makefile
                 src/op/heur/ea/Makefile
                 src/op/exact/Makefile
                 src/op/exact/bac/Makefile
                 ])
AC_OUTPUT

AC_MSG_NOTICE([summary of build options:
               Version:           ${PACKAGE_STRING}
               Debug Build:       ${debug}
               Host type:         ${host}
               Install prefix:    ${prefix}
               Compiler:          ${CC} ${CFLAGS}
               LP solver path:    ${with_cplex}
               Static libraries:  ${enable_static}
               Shared libraries:  ${enable_shared}
               Linker:            ${LD} ${LDFLAGS} ${LIBS}])
