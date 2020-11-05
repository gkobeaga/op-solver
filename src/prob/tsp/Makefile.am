include_HEADERS        = tsp.h
noinst_LTLIBRARIES        = libsolver-prob-tsp.la

libsolver_prob_tsp_la_SOURCES = \
				   env.c  \
				   param.c  \
				   stats.c

libsolver_prob_tsp_la_CPPFLAGS = \
				    -I$(top_srcdir)/src \
				    -I$(top_srcdir)/src/prob

libsolver_prob_tsp_la_LIBADD   =
SUBDIRS  =

SUBDIRS += init
libsolver_prob_tsp_la_LIBADD   += init/libsolver-prob-tsp-init.la
SUBDIRS += heur
libsolver_prob_tsp_la_LIBADD   += heur/libsolver-prob-tsp-heur.la
