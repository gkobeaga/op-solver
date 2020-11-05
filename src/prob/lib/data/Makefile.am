include_HEADERS        = data.h
noinst_LTLIBRARIES        = libsolver-prob-lib-data.la
libsolver_prob_lib_data_la_SOURCES = data.c            \
					io/read.c  \
					nearest/nearest.c  \
					nearest/euclidean.c  \
					nearest/hilbert.c  \
					nearest/banach.c

libsolver_prob_lib_data_la_LIBADD = nearest/kdtree/libsolver-prob-lib-data-kdtree.la

libsolver_prob_lib_data_la_CPPFLAGS = -I$(top_srcdir)/src

SUBDIRS = nearest/kdtree

