#ifndef CP_HEUR_H
#define CP_HEUR_H

#include "../cp.h"

typedef struct cp_heur_param
{
    long time_limit;
    int strat;
#define SOLVER_CP_HEUR_APPR_EA 0
    int improve_sol;
#define SOLVER_CP_ADD_D 0  /* Depending on the distance cost*/
#define SOLVER_CP_ADD_SD 1 /* Depending on distance cost and score increase*/
#define SOLVER_CP_ADD_S 3  /* Depending on the score increase*/
#define SOLVER_CP_ADD_3N 4 /* 3 neigh */
#define SOLVER_CP_ADD_IN 5 /* Intensive */
    int recover_infeas;
#define SOLVER_CP_DROP_D 0  /* Depending on the distance cost*/
#define SOLVER_CP_DROP_SD 1 /* Depending on distance cost and score decrease*/
#define SOLVER_CP_DROP_S 3  /* Dependingon the score decrease*/
    int pop_size;
    int init;
#define SOLVER_CP_INIT_BEST3 0
#define SOLVER_CP_INIT_GREEDY 1
#define SOLVER_CP_INIT_RAND 2
    int select;
#define SOLVER_CP_SEL_BERNOULLI 0 /* select using Bernoully */
    double pinit;                 /* Bernoully p for initial population */
} cp_heur_param;

typedef struct cp_heur_stats
{
    stats_item *total;
} cp_heur_stats;

typedef struct cp_heur_ea_env cp_heur_ea_env;
typedef struct tsp_heur_env tsp_heur_env;
typedef struct cp_heur_env
{
    int verbosity;
    cp_heur_param *param;
    cp_heur_stats *stats;
    cp_heur_ea_env *ea;
    tsp_heur_env *tsp;
    int (*select_best3vertices)(cp_prob *cp, cp_heur_env *env, cp_sol *sol);
    int (*recover_infeas)(cp_prob *cp, cp_heur_env *env, cp_sol *sol);
    int (*local_search)(cp_prob *cp, cp_heur_env *env, cp_sol *sol);
} cp_heur_env;

cp_heur_stats *
cp_create_heur_stats(void);
void
cp_free_heur_stats(cp_heur_stats **stats);

cp_heur_param *
cp_create_heur_param(void);
void
cp_free_heur_param(cp_heur_param **param);

cp_heur_env *
cp_create_heur_env(void);
void
cp_free_heur_env(cp_heur_env **env);
void
cp_conf_heur_env(cp_prob *cp, cp_heur_env *env);

int
cp_improve_heur_in(cp_prob *cp, cp_heur_env *env, cp_sol *sol),
cp_improve_heur_3n(cp_prob *cp, cp_heur_env *env, cp_sol *sol),
cp_recover_heur_infeas(cp_prob *cp, cp_heur_env *env, cp_sol *sol);

int
cp_add_sol_node(cp_prob *cp, cp_heur_env *env, cp_sol *sol, int node),
cp_drop_sol_node(cp_prob *cp, cp_heur_env *env, cp_sol *sol, int node);

int
cp_improve_heur_cycle_length(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol);

int
cp_select_best3nodes(cp_prob *cp, cp_heur_env *env, cp_sol *sol);

#include "ea/ea.h"

#endif
