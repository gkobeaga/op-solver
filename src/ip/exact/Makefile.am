myincludesubdir      = $(pkgincludedir)/ip/
myincludesub_HEADERS = exact.h
noinst_LTLIBRARIES   = libsolver-ip-exact.la

libsolver_ip_exact_la_SOURCES = env.c  \
				param.c  \
				stats.c


libsolver_ip_exact_la_CPPFLAGS = -I$(top_srcdir)/src
SUBDIRS =
libsolver_ip_exact_la_LIBADD   =

SUBDIRS += bac
libsolver_ip_exact_la_LIBADD   += bac/libsolver-ip-exact-bac.la
