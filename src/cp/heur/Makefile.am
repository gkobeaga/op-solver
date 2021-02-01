myincludesubdir      = $(pkgincludedir)/cp/
myincludesub_HEADERS = heur.h
noinst_LTLIBRARIES   = libsolver-cp-heur.la

libsolver_cp_heur_la_SOURCES = \
			       env.c  \
			       param.c  \
			       stats.c  \
			       add.c  \
			       drop.c

libsolver_cp_heur_la_CPPFLAGS = -I$(top_srcdir)/src

libsolver_cp_heur_la_LIBADD   =
SUBDIRS  =

SUBDIRS += ea
libsolver_cp_heur_la_LIBADD   += ea/libsolver-cp-heur-ea.la

