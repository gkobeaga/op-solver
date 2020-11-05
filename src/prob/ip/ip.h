#ifndef IP_H
#define IP_H

#include "../../op-solver.h"
#include "../lp/lp.h"

#define SOLVER_PRICE_MAXPENALTY (0.00001)

typedef struct lp_sol ip_sol;

typedef struct lp_prob lp_prob;
typedef struct ip_prob
{
    solver_graph *graph;
    solver_data *data;
    lp_prob *lp;
    int sense;
    ip_sol *sol;
    double upperboundG;
    double lowerboundG;
    int infeasible;
    void *cuts;
} ip_prob;

typedef struct ip_param
{
    long time_limit;
} ip_param;

typedef struct ip_stats
{
    stats_item *total;
} ip_stats;

typedef struct ip_exact_env ip_exact_env;
typedef struct ip_env
{
    int verbosity;
    ip_param *param;
    ip_stats *stats;
    ip_exact_env *exact;
    void *orig_prob;
    void *orig_env;
} ip_env;

ip_prob *
ip_create_prob(void);
void
ip_free_prob(ip_prob **prob);

ip_env *
ip_create_env(void);
void
ip_free_env(ip_env **env);
ip_param *
ip_create_param(void);
void
ip_free_param(ip_param **param);
ip_stats *
ip_create_stats(void);
void
ip_free_stats(ip_stats **stats);

#define ip_create_sol() lp_data_sol()
#define ip_free_sol(sol) lp_free_sol((sol))
#define ip_update_sol(ip, sol) lp_update_sol((ip)->lp, (sol))
int
cp_is_ipsol_integral_connected(ip_sol *sol);

#include "exact/exact.h"

#endif
