include_HEADERS        = heur.h
noinst_LTLIBRARIES        = libsolver-prob-cp-heur.la

libsolver_prob_cp_heur_la_SOURCES = \
				       env.c  \
				       param.c  \
				       stats.c  \
				       add.c  \
				       drop.c

libsolver_prob_cp_heur_la_CPPFLAGS = \
					-I$(top_srcdir)/src \
					-I$(top_srcdir)/src/prob

libsolver_prob_cp_heur_la_LIBADD   =
SUBDIRS  =

SUBDIRS += ea
libsolver_prob_cp_heur_la_LIBADD   += ea/libsolver-prob-cp-heur-ea.la

