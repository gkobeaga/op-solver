myincludesubdir      = $(pkgincludedir)/graph/
myincludesub_HEADERS = graph.h
lib_LTLIBRARIES      = libsolver-graph.la

libsolver_graph_la_SOURCES = graph.c            \
			     arc_hash.c         \
			     cliques/cliques.c  \
			     cliques/repo.c  \
			     ghtree.c \
			     maxflow.c \
			     mincut.c \
			     mst.c \
			     connect.c

libsolver_graph_la_CPPFLAGS = \
			      -I$(top_srcdir)/src
