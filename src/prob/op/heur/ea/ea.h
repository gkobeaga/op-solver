#ifndef OP_HEUR_EA_H
#define OP_HEUR_EA_H

#include "../heur.h"

typedef struct op_heur_ea_param
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
} op_heur_ea_param;

typedef struct op_heur_ea_stats
{
  stats_item *total;
  stats_item *it;
  stats_item *infeas_recover;
} op_heur_ea_stats;

typedef struct op_heur_ea_env
{
  int verbosity;
  op_heur_ea_param *param;
  op_heur_ea_stats *stats;
} op_heur_ea_env;

op_heur_ea_stats *
op_create_heur_ea_stats(void);
void
op_free_heur_ea_stats(op_heur_ea_stats **stats);

op_heur_ea_param *
op_create_heur_ea_param(void);
void
op_free_heur_ea_param(op_heur_ea_param **param);

op_heur_ea_env *
op_create_heur_ea_env(void);
void
op_free_heur_ea_env(op_heur_ea_env **env);

#define op_opt_heur_ea(prob, env, op, pop,sol) \
  cp_opt_heur_ea((prob), (env), (op), (pop), (sol))
#endif
