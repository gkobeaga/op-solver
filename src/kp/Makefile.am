myincludesubdir      = $(pkgincludedir)/kp/
myincludesub_HEADERS = kp.h
noinst_LTLIBRARIES   = libsolver-kp.la

libsolver_kp_la_SOURCES = prob.c \
			  env.c  \
			  param.c  \
			  stats.c  \
			  sol.c \
			  kp.c \
			  exact/bab.c

libsolver_kp_la_CPPFLAGS = -I$(top_srcdir)/src
