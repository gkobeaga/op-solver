#ifndef TSP_HEUR_H
#define TSP_HEUR_H

#include "../tsp.h"

typedef struct tsp_heur_param
{
    long time_limit;
    int strat;
#define SOLVER_TSP_HEUR_NO 0
#define SOLVER_TSP_HEUR_TWOOPT 1
#define SOLVER_TSP_HEUR_TWOOPT5 2
#define SOLVER_TSP_HEUR_THREEOPT 3
#define SOLVER_TSP_HEUR_LINKERN 4
    int two_and_a_half_opt;
} tsp_heur_param;

typedef struct stats_item stats_item;
typedef struct tsp_heur_stats
{
    stats_item *total;
} tsp_heur_stats;

typedef struct tsp_heur_linkern_env tsp_heur_linkern_env;
typedef struct tsp_heur_env
{
    int verbosity;
    tsp_heur_param *param;
    tsp_heur_stats *stats;
    tsp_heur_linkern_env *linkern;
    int node_a;
    int node_b;
    int node_c;
    int ab_dist;
    int ac_dist;
    int cd_dist;
} tsp_heur_env;

tsp_heur_stats *
tsp_create_heur_stats(void);
void
tsp_free_heur_stats(tsp_heur_stats **stats);

tsp_heur_param *
tsp_create_heur_param(void);
void
tsp_free_heur_param(tsp_heur_param **param);

tsp_heur_env *
tsp_create_heur_env(void);
void
tsp_free_heur_env(tsp_heur_env **env);

int
tsp_opt_heur(tsp_prob *tsp, tsp_heur_env *env, tsp_sol *sol),
tsp_opt_heur_threeopt(tsp_prob *tsp, tsp_heur_env *env, tsp_sol *sol),
tsp_opt_heur_twoopt(tsp_prob *tsp, tsp_heur_env *heur_env, tsp_sol *sol),
tsp_opt_heur_twoopt_banach(tsp_prob *tsp, tsp_heur_env *heur_env, tsp_sol *sol);

#include "linkern/linkern.h"

#endif
