#ifndef TSP_INIT_H
#define TSP_INIT_H

typedef struct tsp_init_param
{
    long time_limit;
    int init;
#define TSP_INIT_RANDOM 0
#define TSP_INIT_NEIGHBOR 1
#define TSP_INIT_GREEDY 2
#define TSP_INIT_BORUVKA 3
#define TSP_INIT_QBORUVKA 4
} tsp_init_param;

typedef struct stats_item stats_item;
typedef struct tsp_init_stats
{
    stats_item *total;
} tsp_init_stats;

typedef struct tsp_init_env
{
    int verbosity;
    tsp_init_param *param;
    tsp_init_stats *stats;
} tsp_init_env;

tsp_init_stats *
tsp_create_init_stats(void);
void
tsp_free_init_stats(tsp_init_stats **stats);

tsp_init_param *
tsp_create_init_param(void);
void
tsp_free_init_param(tsp_init_param **param);

tsp_init_env *
tsp_create_init_env(void);
void
tsp_free_init_env(tsp_init_env **env);

typedef struct cp_prob tsp_prob;
typedef struct cp_sol tsp_sol;
typedef struct cp_pop tsp_pop;

int
tsp_init_sol(tsp_prob *tsp, tsp_init_env *env, tsp_sol *sol);

int
tsp_init_pop(tsp_prob *tsp, tsp_init_env *env, tsp_pop *pop);

int
tsp_init_sol_random(tsp_prob *tsp, tsp_init_env *env, tsp_sol *sol);
#endif
