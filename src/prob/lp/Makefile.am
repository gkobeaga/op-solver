include_HEADERS        = lp.h
noinst_LTLIBRARIES        = libsolver-prob-lp.la

libsolver_prob_lp_la_LIBADD = \
				 $(CPLEX_PATH)/lib/x86-64_sles10_4.1/static_pic/libcplex.a \
				 -ldl -lm -lpthread
libsolver_prob_lp_la_SOURCES = prob.c            \
				  data.c         \
				  env.c  \
				  param.c  \
				  stats.c  \
				  sol.c \
				  lib/cplex.c

libsolver_prob_lp_la_CPPFLAGS = \
				   -I$(top_srcdir)/src \
				   -I$(top_srcdir)/src/prob \
				   -I$(CPLEX_PATH)/include
