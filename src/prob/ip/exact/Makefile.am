include_HEADERS        = exact.h
noinst_LTLIBRARIES        = libsolver-prob-ip-exact.la

libsolver_prob_ip_exact_la_SOURCES = env.c  \
					param.c  \
					stats.c


libsolver_prob_ip_exact_la_CPPFLAGS = \
					 -I$(top_srcdir)/src \
					 -I$(top_srcdir)/src/prob

SUBDIRS =
libsolver_prob_ip_exact_la_LIBADD   =

SUBDIRS += bac
libsolver_prob_ip_exact_la_LIBADD   += bac/libsolver-prob-ip-exact-bac.la
