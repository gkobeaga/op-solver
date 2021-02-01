#ifndef KP_H
#define KP_H

#include "../op-solver.h"

typedef struct kp_sol
{
    int tot_n;
    int *selected;
    int *sposition;
    double val;
    double weight;
    int ns;
} kp_sol;

typedef struct kp_prob
{
    int n;
    int nne; /* Number of non-excluded items (price_i < cap) */
    double *w;
    double *p;
    double cap;
    double tot_profit;
    int sol_stat;
    int *ord_ind;
    double *ord_w;
    double *ord_p;
    double upperbound;
    double *min_w;
    kp_sol *sol;
} kp_prob;

typedef struct kp_param
{
    long time_limit;
    int exact_tech;
    int check_input;
    int reorder_items;
    double epsilon;
    double p_epsilon;
    double w_epsilon;
    const char *stats_file;
} kp_param;

typedef struct kp_stats
{
    struct stats_item *total;
    stats_item *init_sol;
    stats_item *exact;
    stats_item *branch_node;
} kp_stats;

typedef struct kp_env
{
    int verbosity;
    kp_param *param;
    kp_stats *stats;
    double branch_cap;
    double branch_profit;
    int p_integer;
    int w_integer;
    int counter;
    int max_counter;
    int j;
    int r;
#define KP_EXACT_NONE 0
#define KP_EXACT_BAB 1
} kp_env;

kp_prob *
kp_create_prob(solver_data *data);
void
kp_erase_prob(kp_prob *kp);
void
kp_free_prob(kp_prob **kp);
int
kp_read_prob(kp_prob *kp, kp_env *env, const char *fname, int flags);

int
kp_opt(kp_prob *kp, kp_env *env, kp_sol *sol);

kp_env *
kp_create_env(void);
void
kp_free_env(kp_env **env);
kp_param *
kp_create_param(void);
void
kp_free_param(kp_param **param);
kp_stats *
kp_create_stats(void);
void
kp_free_stats(kp_stats **stats);

int
kp_get_mincover(kp_prob *kp, kp_sol *sol);

kp_sol *
kp_create_sol(kp_prob *kp);
void
kp_copy_sol(kp_sol *insol, kp_sol *outsol);
void
kp_erase_sol(kp_sol *sol);
void
kp_free_sol(kp_sol **sol);

int
kp_start_sol(kp_prob *kp, kp_env *env, kp_sol *sol);
int
kp_opt_bab(kp_prob *kp, kp_env *env, kp_sol *sol);
#endif
