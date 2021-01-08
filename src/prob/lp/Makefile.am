include_HEADERS    = lp.h
noinst_LTLIBRARIES = libsolver-prob-lp.la

libsolver_prob_lp_la_LIBADD =  @LP_SOLVER_LIB@ -ldl -lm -lpthread
libsolver_prob_lp_la_SOURCES = prob.c            \
			       data.c         \
			       env.c  \
			       param.c  \
			       stats.c  \
			       sol.c
if HAVE_LP_SOLVER
libsolver_prob_lp_la_SOURCES+=lib/cplex.c
endif

libsolver_prob_lp_la_CPPFLAGS = \
				-I$(top_srcdir)/src \
				-I$(top_srcdir)/src/prob \
				-I@LP_SOLVER_INCLUDE@
