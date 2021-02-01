myincludesubdir      = $(pkgincludedir)/tsp/
myincludesub_HEADERS = heur.h
noinst_LTLIBRARIES   = libsolver-tsp-heur.la

#myincludesubdir      = $(pkgincludedir)/tsp/heur/
#myincludesub_HEADERS = linkern/linkern.h

libsolver_tsp_heur_la_SOURCES = \
					env.c  \
					param.c  \
					stats.c  \
					opt.c  \
					kopt/ls.c \
					linkern/stats.c  \
					linkern/env.c  \
					linkern/param.c  \
					linkern/flip_two.c  \
					linkern/linkern.c

libsolver_tsp_heur_la_CPPFLAGS = -I$(top_srcdir)/src

