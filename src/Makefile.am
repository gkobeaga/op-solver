pkginclude_HEADERS= op-solver.h

SUBDIRS = lib data graph
SUBDIRS += cp ip kp lp op tsp

lib_LTLIBRARIES     = libopsolver.la
libopsolver_la_LIBADD = \
			data/libsolver-data.la \
			graph/libsolver-graph.la \
			lib/libsolver-lib.la \
			cp/libsolver-cp.la \
			kp/libsolver-kp.la \
			tsp/libsolver-tsp.la

if HAVE_LP_SOLVER
libopsolver_la_LIBADD += \
			 ip/libsolver-ip.la \
			 lp/libsolver-lp.la \
			 @LP_SOLVER_LIB@ -lpthread -lgmp
endif
libopsolver_la_SOURCES =


bin_PROGRAMS = op-solver
op_solver_SOURCES = op-solver.c cmd/opt.c
op_solver_LDADD = libopsolver.la -lm
op_solver_CPPFLAGS = -I$(top_srcdir)/src
