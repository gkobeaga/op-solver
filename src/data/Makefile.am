include_HEADERS        = data.h
noinst_LTLIBRARIES        = libsolver-data.la
libsolver_data_la_SOURCES = data.c            \
			    io/read.c  \
			    nearest/nearest.c  \
			    nearest/euclidean.c  \
			    nearest/hilbert.c  \
			    nearest/banach.c

libsolver_data_la_LIBADD = nearest/kdtree/libsolver-data-kdtree.la

libsolver_data_la_CPPFLAGS = -I$(top_srcdir)/src

SUBDIRS = nearest/kdtree

