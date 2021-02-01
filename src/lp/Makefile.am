myincludesubdir      = $(pkgincludedir)/lp/
myincludesub_HEADERS = lp.h
noinst_LTLIBRARIES   = libsolver-lp.la

libsolver_lp_la_LIBADD =  @LP_SOLVER_LIB@ -ldl -lm -lpthread
libsolver_lp_la_SOURCES = prob.c \
			  data.c \
			  env.c \
			  param.c \
			  stats.c \
			  sol.c
if HAVE_LP_SOLVER
if HAVE_LP_SOLVER_CPLEX
libsolver_lp_la_SOURCES+=lib/cplex.c
endif
if HAVE_LP_SOLVER_DUMMY
libsolver_lp_la_SOURCES+=lib/dummy.c
endif
endif

libsolver_lp_la_CPPFLAGS = -I$(top_srcdir)/src \
			   -I@LP_SOLVER_INCLUDE@
