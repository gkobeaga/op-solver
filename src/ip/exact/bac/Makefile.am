myincludesubdir      = $(pkgincludedir)/cp/exact
myincludesub_HEADERS = bac.h
noinst_LTLIBRARIES   = libsolver-ip-exact-bac.la

libsolver_ip_exact_bac_la_SOURCES = env.c  \
					    param.c  \
					    stats.c  \
					    branch/branch.c \
					    branch/history.c \
					    branch/dfs.c

libsolver_ip_exact_bac_la_CPPFLAGS = -I$(top_srcdir)/src
