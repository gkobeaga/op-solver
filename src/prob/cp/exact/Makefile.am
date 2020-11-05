include_HEADERS        = exact.h
lib_LTLIBRARIES        = libsolver-prob-cp-exact.la

libsolver_prob_cp_exact_la_SOURCES = \
					env.c  \
					param.c  \
					stats.c  \
					opt.c


libsolver_prob_cp_exact_la_CPPFLAGS = \
					 -I$(top_srcdir)/src \
					 -I$(top_srcdir)/src/prob

libsolver_prob_cp_exact_la_LIBADD   =
SUBDIRS  =

SUBDIRS += bac
libsolver_prob_cp_exact_la_LIBADD   += bac/libsolver-prob-cp-exact-bac.la
