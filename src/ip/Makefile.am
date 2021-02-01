myincludesubdir      = $(pkgincludedir)/ip/
myincludesub_HEADERS = ip.h
noinst_LTLIBRARIES   = libsolver-ip.la

libsolver_ip_la_SOURCES = prob.c \
			  env.c  \
			  param.c  \
			  stats.c  \
			  sol.c


libsolver_ip_la_CPPFLAGS = -I$(top_srcdir)/src
SUBDIRS =
libsolver_ip_la_LIBADD   =

SUBDIRS += exact
libsolver_ip_la_LIBADD   += exact/libsolver-ip-exact.la
