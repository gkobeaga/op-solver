include_HEADERS        = bac.h
lib_LTLIBRARIES        = libsolver-prob-cp-exact-bac.la

libsolver_prob_cp_exact_bac_la_SOURCES = \
					 env.c  \
					 param.c  \
					 stats.c  \
					 opt.c  \
					 init.c  \
					 dual.c  \
					 pricing.c  \
					 xheur.c  \
					 branch/simple.c  \
					 lib/lp.c  \
					 cut/age.c  \
					 cut/cut.c  \
					 cut/hash.c  \
					 cut/loop.c  \
					 cut/skeleton.c  \
					 cut/logical.c  \
					 cut/comb/blossom_fast.c  \
					 cut/comb/blossom_ghfast.c  \
					 cut/comb/blossom_mst.c  \
					 cut/cover/cycle.c  \
					 cut/cover/edge.c  \
					 cut/cover/vertex.c  \
					 cut/path/path.c  \
					 cut/sec/sec.c  \
					 cut/sec/hong.c  \
					 cut/sec/gomoryhu.c \
					 shrinking/shrinking.c  \
					 shrinking/c1.c  \
					 shrinking/c1c2.c  \
					 shrinking/c1c2c3.c  \
					 shrinking/s1.c

libsolver_prob_cp_exact_bac_la_CPPFLAGS = \
					  -I$(top_srcdir)/src \
					  -I$(top_srcdir)/src/prob
