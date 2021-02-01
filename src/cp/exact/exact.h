#ifndef CP_EXACT_H
#define CP_EXACT_H

#include "../cp.h"

typedef struct cp_exact_param
{
    long time_limit;
    int appr;
} cp_exact_param;

typedef struct cp_exact_stats
{
    stats_item *total;
} cp_exact_stats;

typedef struct cp_exact_bac_env cp_exact_bac_env;
typedef struct cp_exact_env
{
    int verbosity;
    cp_exact_param *param;
    cp_exact_stats *stats;
    cp_exact_bac_env *bac;
} cp_exact_env;

cp_exact_stats *
cp_create_exact_stats(void);
void
cp_free_exact_stats(cp_exact_stats **stats);

cp_exact_param *
cp_create_exact_param(void);
void
cp_free_exact_param(cp_exact_param **param);

cp_exact_env *
cp_create_exact_env(void);
void
cp_free_exact_env(cp_exact_env **env);
int
cp_parse_exact_args(int argc, char *argv[], cp_exact_env *env);

int
cp_opt_exact(cp_prob *cp, cp_exact_env *exact_env, cp_sol *sol);

#include "bac/bac.h"

#endif
