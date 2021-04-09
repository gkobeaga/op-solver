myincludesubdir      = $(pkgincludedir)/tsp/
myincludesub_HEADERS = heur.h
noinst_LTLIBRARIES   = libsolver-tsp-heur.la

libsolver_tsp_heur_la_SOURCES = \
				env.c  \
				param.c  \
				stats.c  \
				opt.c  \
				kopt/ls.c

libsolver_tsp_heur_la_CPPFLAGS = -I$(top_srcdir)/src

libsolver_tsp_heur_la_LIBADD   =
SUBDIRS  =

SUBDIRS += linkern
libsolver_tsp_heur_la_LIBADD   += linkern/libsolver-tsp-heur-linkern.la

