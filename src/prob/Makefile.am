noinst_LTLIBRARIES		= libsolver-prob.la
libsolver_prob_la_SOURCES =
libsolver_prob_la_LIBADD   =
SUBDIRS =

SUBDIRS += lib
libsolver_prob_la_LIBADD   += lib/libsolver-prob-lib.la

if HAVE_LP_SOLVER
# Linear Problem
SUBDIRS += lp
libsolver_prob_la_LIBADD   += lp/libsolver-prob-lp.la
endif

if HAVE_LP_SOLVER
# Integer Problem
SUBDIRS += ip
libsolver_prob_la_LIBADD   += ip/libsolver-prob-ip.la
endif

# Knapsack Problem
SUBDIRS += kp
libsolver_prob_la_LIBADD   += kp/libsolver-prob-kp.la

# Cycle Problem
SUBDIRS += cp
libsolver_prob_la_LIBADD   += cp/libsolver-prob-cp.la

# Travelling Salesperson Problem
SUBDIRS += tsp
libsolver_prob_la_LIBADD   += tsp/libsolver-prob-tsp.la

#Orienteering Problem
SUBDIRS += op
