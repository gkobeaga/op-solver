include_HEADERS        = ip.h
noinst_LTLIBRARIES     = libsolver-prob-ip.la

libsolver_prob_ip_la_SOURCES = prob.c            \
				  env.c  \
				  param.c  \
				  stats.c  \
				  sol.c


libsolver_prob_ip_la_CPPFLAGS = \
				   -I$(top_srcdir)/src \
				   -I$(top_srcdir)/src/prob

SUBDIRS =
libsolver_prob_ip_la_LIBADD   =

SUBDIRS += exact
libsolver_prob_ip_la_LIBADD   += exact/libsolver-prob-ip-exact.la
