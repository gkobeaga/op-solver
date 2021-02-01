#ifndef OP_INIT_H
#define OP_INIT_H

#include "../op.h"

typedef struct op_init_param
{
    long time_limit;
    int init;
#define SOLVER_OP_INIT_BEST3 0
#define SOLVER_OP_INIT_GREEDY 1
#define SOLVER_OP_INIT_RAND 2
    int select;
#define SOLVER_OP_SEL_BERNOULLI 0
    double pinit;
} op_init_param;

typedef struct op_init_stats
{
    stats_item *total;
} op_init_stats;

typedef struct op_init_env
{
    int verbosity;
    op_init_param *param;
    op_init_stats *stats;
} op_init_env;

op_init_stats *
op_create_init_stats(void);
void
op_free_init_stats(op_init_stats **stats);

op_init_param *
op_create_init_param(void);
void
op_free_init_param(op_init_param **param);

op_init_env *
op_create_init_env(void);
void
op_free_init_env(op_init_env **env);

#define op_init_sol(prob, env, op, sol) cp_init_sol((prob), (env), (op), (sol))
#define op_init_pop(prob, env, op, pop) cp_init_pop((prob), (env), (op), (pop))
#endif
