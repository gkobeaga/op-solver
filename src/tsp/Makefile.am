myincludesubdir      = $(pkgincludedir)/tsp/
myincludesub_HEADERS = tsp.h
noinst_LTLIBRARIES   = libsolver-tsp.la

libsolver_tsp_la_SOURCES = \
			   env.c  \
			   param.c  \
			   stats.c

libsolver_tsp_la_CPPFLAGS = -I$(top_srcdir)/src

libsolver_tsp_la_LIBADD   =
SUBDIRS  =

SUBDIRS += init
libsolver_tsp_la_LIBADD   += init/libsolver-tsp-init.la
SUBDIRS += heur
libsolver_tsp_la_LIBADD   += heur/libsolver-tsp-heur.la
