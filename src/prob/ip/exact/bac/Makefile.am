include_HEADERS        = bac.h
noinst_LTLIBRARIES        = libsolver-prob-ip-exact-bac.la

libsolver_prob_ip_exact_bac_la_SOURCES = env.c  \
					    param.c  \
					    stats.c  \
					    branch/branch.c \
					    branch/history.c \
					    branch/dfs.c

libsolver_prob_ip_exact_bac_la_CPPFLAGS = \
					     -I$(top_srcdir)/src \
					     -I$(top_srcdir)/src/prob
