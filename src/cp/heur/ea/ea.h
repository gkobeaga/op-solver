#ifndef CP_HEUR_EA_H
#define CP_HEUR_EA_H

#include "../heur.h"

typedef struct cp_heur_ea_param
{
    long time_limit;
    int it_lim;
    int pop_size;
    int len_improve1;
    int len_improve2;
    double pmut;
    int nparsel;
    double pinit;
    int d2d;
    int pop_stop;
} cp_heur_ea_param;

typedef struct cp_heur_ea_stats
{
    stats_item *total;
    stats_item *it;
    stats_item *infeas_recover;
} cp_heur_ea_stats;

typedef struct cp_heur_ea_env
{
    int verbosity;
    cp_heur_ea_param *param;
    cp_heur_ea_stats *stats;
    cp_heur_env *heur;
} cp_heur_ea_env;

cp_heur_ea_stats *
cp_create_heur_ea_stats(void);
void
cp_free_heur_ea_stats(cp_heur_ea_stats **stats);

cp_heur_ea_param *
cp_create_heur_ea_param(void);
void
cp_free_heur_ea_param(cp_heur_ea_param **param);

cp_heur_ea_env *
cp_create_heur_ea_env(void);
void
cp_free_heur_ea_env(cp_heur_ea_env **env);

int
cp_opt_heur_ea(cp_prob *cp, cp_heur_env *ea_env, cp_pop *pop, cp_sol *sol);
int
cp_heur_ea_selection(cp_prob *cp, cp_heur_ea_env *heur_ea_env, cp_pop *pop,
                     int *parents);
int
cp_heur_ea_crossover(cp_prob *cp, cp_heur_ea_env *heur_env, cp_sol *par0,
                     cp_sol *par1, cp_sol *child);
int
cp_heur_ea_mutation(cp_prob *cp, cp_heur_ea_env *heur_env, cp_sol *sol);
#endif
