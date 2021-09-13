#ifndef SOLVER_H
#define SOLVER_H

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <gmp.h>

#include <inttypes.h>

#include "data/data.h"
#include "graph/graph.h"

#ifdef __cplusplus
extern "C"
{
#endif

    /* library version numbers: */
#define SOLVER_MAJOR_VERSION 0
#define SOLVER_MINOR_VERSION 1

#define SOLVER_UNDEF 0  /* integer solution is undefined */
#define SOLVER_OPT 1    /* integer solution is optimal */
#define SOLVER_FEAS 2   /* integer solution is feasible */
#define SOLVER_NOFEAS 3 /* no integer solution exists */

#define SOLVER_SOLVE_TSP 0
#define SOLVER_SOLVE_OP 1

#define SOLVER_OPT_SENSE_MIN 0
#define SOLVER_OPT_SENSE_MAX 1

#define HASH_UPDATE_NAME 0
#define HASH_UPDATE_OPD0 1
#define HASH_UPDATE_OPSCORE 2

#define SOLVER_ONEMINUS 0.9999999
#define SOLVER_ZEROPLUS 0.0000001

#define SOLVER_MAXDOUBLE 1.0e30

#define SOLVER_VERBOSITY_OFF 0
#define SOLVER_VERBOSITY_ERROR 1
#define SOLVER_VERBOSITY_INFO 2
#define SOLVER_VERBOSITY_DEBUG 3

    extern unsigned long __seed__;

#ifndef SOLVER_MAXINT
#define SOLVER_MAXINT ((int)(~(((unsigned)1) << ((8 * sizeof(int)) - 1))))
#endif

#define SWAP(a, b, t) (((t) = (a)), ((a) = (b)), ((b) = (t)))

#define SOLVER_CTERM_RESET "\033[0m"
#define SOLVER_CTERM_BOLDGREEN "\033[1m\033[32m"
#define SOLVER_CTERM_BOLDRED "\033[1m\033[31m"
#define SOLVER_CTERM_BOLDCYAN "\033[1m\033[36m"

#define __FILENAME__ (strrchr(__FILE__, 's') + 4)

#define check_rval(rval, msg, label)                                           \
    {                                                                          \
        if ((rval))                                                            \
        {                                                                      \
            fprintf(stderr,                                                    \
                    SOLVER_CTERM_BOLDRED "[%s:%d] " SOLVER_CTERM_BOLDCYAN      \
                                         "%s: " SOLVER_CTERM_RESET             \
                                         "rval %d - %s\n",                     \
                    __FILENAME__, __LINE__, __FUNCTION__, rval, (msg));        \
            goto label;                                                        \
        }                                                                      \
    }

#define check_assert(eval, msg, label)                                         \
    {                                                                          \
        if (!(eval))                                                           \
        {                                                                      \
            fprintf(stderr,                                                    \
                    SOLVER_CTERM_BOLDRED "[%s:%d] " SOLVER_CTERM_BOLDCYAN      \
                                         "%s: " SOLVER_CTERM_RESET             \
                                         "Assertion `%s` failed. - %s\n",      \
                    __FILENAME__, __LINE__, __FUNCTION__, #eval, (msg));       \
            rval = 1;                                                          \
            goto label;                                                        \
        }                                                                      \
        else                                                                   \
            rval = 0;                                                          \
    }

#define check_null(item, msg, label)                                           \
    {                                                                          \
        if ((!item))                                                           \
        {                                                                      \
            fprintf(stderr,                                                    \
                    SOLVER_CTERM_BOLDRED "[%s:%d] " SOLVER_CTERM_BOLDCYAN      \
                                         "%s: " SOLVER_CTERM_RESET             \
                                         "%s is NULL - %s \n",                 \
                    __FILENAME__, __LINE__, __FUNCTION__, (#item), (msg));     \
            rval = 1;                                                          \
            goto label;                                                        \
        }                                                                      \
    }

#define realloc_scale(ptr, pnnum, count, scale)                                \
    {                                                                          \
        __typeof__(ptr) p;                                                     \
        int newsize = (int)(((double)pnnum) * (scale));                        \
        newsize     = (newsize < pnnum + 1000 ? (pnnum + 1000) : newsize);     \
        newsize     = (newsize < count ? (count) : newsize);                   \
        p =                                                                    \
        ((__typeof__(ptr))realloc(ptr, newsize * sizeof(__typeof__(ptr))));    \
        pnnum = ((p) ? newsize : 0);                                           \
        ptr   = p;                                                             \
    }

    int __solver_opt(int argc, char *argv[]);

    void sort_partial(int *arr, int l, int r, int m, double *coord);

    typedef struct stats_item
    {
        char name[40];
        struct timespec current_start;
        struct timespec current_end;
        struct timespec last_start;
        struct timespec last_end;
        long last_time;
        long total_time;
        int count_total;
        int count_success;
        int count_active;
        int last_success;
        int last_total;
        int active;
    } stats_item;

    stats_item *stats_create(const char *name);

    int stats_start(stats_item *t), stats_stop(stats_item *t, int count),
    stats_stop_if_active(stats_item *t, int count),
    stats_suspend(stats_item *t), stats_resume(stats_item *t),
    stats_print(stats_item *t);

    long stats_get_current_time(stats_item *t),
    stats_get_total_time(stats_item *t);

    /* Prime Numbers */
    int prime_check(unsigned int p);
    unsigned int prime_next(unsigned int x);

    typedef struct solver_dheap
    {
        double *key;
        int *entry;
        int *loc;
        int total_space;
        int size;
    } solver_dheap;

    /* dheap*/
    void dheap_free(solver_dheap **h), dheap_delete(solver_dheap *h, int i),
    dheap_changekey(solver_dheap *h, int i, double newkey);

    solver_dheap *dheap_create(int k);
    int dheap_resize(solver_dheap *h, int newsize),
    dheap_findmin(solver_dheap *h), dheap_deletemin(solver_dheap *h),
    dheap_insert(solver_dheap *h, int i);

    /* Random */
    int rng_bernoulli(double p);
    void rng_choose(void *candidates, int n, void *selected, int k,
                    size_t size);
    void rng_choose_p(void *candidates, double *p, int ncand, void *selected,
                      int nsel, size_t size);

#include "lp/lp.h"
#include "ip/ip.h"
#include "kp/kp.h"
#include "cp/cp.h"
#include "tsp/tsp.h"
#include "op/op.h"

#ifdef __cplusplus
}
#endif

#endif
