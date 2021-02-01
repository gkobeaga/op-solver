#ifndef TSP_H
#define TSP_H

#include "../cp/cp.h"

typedef struct cp_sol tsp_sol;
typedef struct cp_pop tsp_pop;
typedef struct cp_prob tsp_prob;

typedef struct tsp_param
{
    long time_limit;
    int check_input;
    int reorder_items;
    int exact; /* Find exact solution */
#define SOLVER_TSP_EXACT_NONE 0
#define SOLVER_TSP_EXACT_BAC 1
    int heur;
#define SOLVER_TSP_HEUR_NONE 0
#define SOLVER_TSP_HEUR_EA 1
} tsp_param;

typedef struct stats_item stats_item;
typedef struct tsp_stats
{
    stats_item *total;
} tsp_stats;

typedef struct tsp_init_env tsp_init_env;
typedef struct tsp_heur_env tsp_heur_env;
typedef struct tsp_exact_env tsp_exact_env;
typedef struct tsp_env
{
    int verbosity;
    tsp_param *param;
    tsp_stats *stats;
    tsp_init_env *init;
    tsp_heur_env *heur;
    tsp_exact_env *exact;
} tsp_env;

tsp_env *
tsp_create_env(void);
void
tsp_free_env(tsp_env **env);
tsp_param *
tsp_create_param(void);
void
tsp_free_param(tsp_param **param);
tsp_stats *
tsp_create_stats(void);
void
tsp_free_stats(tsp_stats **stats);

#define tsp_create_prob(data) cp_create_prob((data))
#define tsp_erase_prob(prob) cp_erase_prob((prob))
#define tsp_free_prob(prob) cp_free_prob((prob))

int
tsp_opt(tsp_prob *tsp, tsp_env *env, tsp_sol *sol);

#define tsp_create_sol(prob) cp_create_sol((prob))
#define tsp_copy_sol(insol, outsol) cp_copy_sol((insol), (outsol))
#define tsp_erase_sol(sol) cp_erase_sol((sol))
#define tsp_free_sol(sol) cp_free_sol((sol))

#define tsp_create_pop(prob, size) cp_create_pop((prob), (size))
#define tsp_set_pop_sol(pop, sol, pos) cp_set_pop_sol((pop), (sol), (pos))
#define tsp_update_pop(pop) cp_update_pop((pop))
#define tsp_erase_pop(pop) cp_erase_pop((pop))
#define tsp_free_pop(pop) cp_free_pop((pop))

#endif
