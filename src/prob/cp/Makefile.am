include_HEADERS        = cp.h
noinst_LTLIBRARIES     = libsolver-prob-cp.la


libsolver_prob_cp_la_SOURCES = prob.c            \
			       env.c  \
			       param.c  \
			       stats.c  \
			       sol.c \
			       pop.c \
			       opt.c \
			       lib/tsp.c

libsolver_prob_cp_la_CPPFLAGS = \
				-I$(top_srcdir)/src \
				-I$(top_srcdir)/src/prob

libsolver_prob_cp_la_LIBADD   =
SUBDIRS  =

SUBDIRS += init
libsolver_prob_cp_la_LIBADD   += init/libsolver-prob-cp-init.la
SUBDIRS += heur
libsolver_prob_cp_la_LIBADD   += heur/libsolver-prob-cp-heur.la
if HAVE_LP_SOLVER
SUBDIRS += exact
libsolver_prob_cp_la_CFLAGS   = -DHAVE_LP_SOLVER
libsolver_prob_cp_la_LIBADD   += exact/libsolver-prob-cp-exact.la
endif
