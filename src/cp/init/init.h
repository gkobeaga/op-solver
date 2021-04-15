#ifndef CP_INIT_H
#define CP_INIT_H

#include "../cp.h"

typedef struct cp_init_param
{
    long time_limit;
    int init;
#define SOLVER_CP_INIT_BEST3 0
#define SOLVER_CP_INIT_GREEDY 1
#define SOLVER_CP_INIT_RAND 2
    int select;
#define SOLVER_CP_SEL_BERNOULLI 0 /* select using Bernoully */
    double pinit;                 /* Bernoully p for initial population */
} cp_init_param;

typedef struct cp_init_stats
{
    stats_item *total;
    int write;
    const char *file;
} cp_init_stats;

typedef struct cp_init_env
{
    int verbosity;
    cp_init_param *param;
    cp_init_stats *stats;
} cp_init_env;

cp_init_stats *
cp_create_init_stats(void);
void
cp_free_init_stats(cp_init_stats **stats);
int
cp_write_init_stats(cp_prob *cp, cp_init_env *env);

cp_init_param *
cp_create_init_param(void);
void
cp_free_init_param(cp_init_param **param);

cp_init_env *
cp_create_init_env(void);
void
cp_free_init_env(cp_init_env **env);

int
cp_select_sol_vertices(cp_prob *cp, cp_heur_env *env, cp_sol *sol);

int
cp_init_sol(cp_prob *cp, cp_heur_env *env, cp_sol *sol);

int
cp_init_pop(cp_prob *cp, cp_heur_env *env, cp_pop *pop);

#endif
