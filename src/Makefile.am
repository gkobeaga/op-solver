include_HEADERS		= op-solver.h

SUBDIRS = lib prob

lib_LTLIBRARIES		= libsolver.la
libsolver_la_LIBADD  = lib/libsolver-lib.la \
			  prob/libsolver-prob.la \
			  $(CPLEX_PATH)/lib/x86-64_sles10_4.1/static_pic/libcplex.a -lpthread -lgmp
libsolver_la_SOURCES =


bin_PROGRAMS = op-solver
op_solver_SOURCES = op-solver.c cmd/opt.c
op_solver_LDADD = libsolver.la -lm
op_solver_CPPFLAGS = -I$(top_srcdir)/src
