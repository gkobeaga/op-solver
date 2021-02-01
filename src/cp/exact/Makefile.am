myincludesubdir      = $(pkgincludedir)/cp/
myincludesub_HEADERS = exact.h
lib_LTLIBRARIES      = libsolver-cp-exact.la

libsolver_cp_exact_la_SOURCES = \
				env.c  \
				param.c  \
				stats.c  \
				opt.c


libsolver_cp_exact_la_CPPFLAGS = -I$(top_srcdir)/src

libsolver_cp_exact_la_LIBADD   =
SUBDIRS  =

SUBDIRS += bac
libsolver_cp_exact_la_LIBADD   += bac/libsolver-cp-exact-bac.la
