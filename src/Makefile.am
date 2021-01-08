include_HEADERS		= op-solver.h

SUBDIRS = lib data graph prob

lib_LTLIBRARIES     = libsolver.la
libsolver_la_LIBADD = data/libsolver-data.la \
		      graph/libsolver-graph.la \
		      lib/libsolver-lib.la \
		      prob/libsolver-prob.la

if HAVE_LP_SOLVER
libsolver_la_LIBADD += @LP_SOLVER_LIB@ -lpthread -lgmp
endif
libsolver_la_SOURCES =


bin_PROGRAMS = op-solver
op_solver_SOURCES = op-solver.c cmd/opt.c
op_solver_LDADD = libsolver.la -lm
op_solver_CPPFLAGS = -I$(top_srcdir)/src
