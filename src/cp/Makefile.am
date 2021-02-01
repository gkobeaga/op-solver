myincludesubdir      = $(pkgincludedir)/cp/
myincludesub_HEADERS = cp.h
noinst_LTLIBRARIES   = libsolver-cp.la


libsolver_cp_la_SOURCES = prob.c \
			  env.c  \
			  param.c  \
			  stats.c  \
			  sol.c \
			  pop.c \
			  opt.c \
			  lib/tsp.c

libsolver_cp_la_CPPFLAGS = -I$(top_srcdir)/src

libsolver_cp_la_LIBADD   =
SUBDIRS  =

SUBDIRS += init
libsolver_cp_la_LIBADD   += init/libsolver-cp-init.la
SUBDIRS += heur
libsolver_cp_la_LIBADD   += heur/libsolver-cp-heur.la
if HAVE_LP_SOLVER
SUBDIRS += exact
libsolver_cp_la_CFLAGS   = -DHAVE_LP_SOLVER
libsolver_cp_la_LIBADD   += exact/libsolver-cp-exact.la
endif
