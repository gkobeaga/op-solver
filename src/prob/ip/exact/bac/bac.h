#ifndef IP_EXACT_BAC_H
#define IP_EXACT_BAC_H

#include "../../../ip/ip.h"
#include "../exact.h"

#define SOLVER_IP_MAX_COLS (100000)
#define SOLVER_IP_PRICE_MAX_ADD (200)

#define SOLVER_IP_BAC_PRICE_MAXPENALTY (0.00001)
#define SOLVER_IP_BAC_PHASE1_MAXPENALTY (0.001)

#define SOLVER_IP_BAC_PRICE_RCTHRESH (0.00001)
#define SOLVER_IP_BAC_PHASE1_RCTHRESH (0.001)

#define SOLVER_IP_BAC_PRICE_GEN (20000)
#define SOLVER_IP_BAC_PRICE_GEN_FACTOR (3)
#define SOLVER_IP_BAC_PRICE_POOL (1000)

#define SOLVER_IP_BAC_STORE_BATCH (50)
#define SOLVER_IP_BAC_CUT_BATCH (250)

#define SOLVER_IP_BAC_DUAL_DUST (0.001)
#define SOLVER_IP_BAC_SLACK_DUST (0.001)
#define SOLVER_IP_BAC_VAR_DUST (0.001)

#define SOLVER_IP_BAC_MAX_CUT_AGE (5)
#define SOLVER_IP_BAC_MAX_VAR_AGE (100)

#define SOLVER_IP_BAC_MIN_VIOL (0.000001)

typedef struct ip_exact_bac_param
{
    long time_limit;
    double nexttol;
    int branch_strat;
#define SOLVER_IP_SEARCH_DFS 0
    int branch_select;
#define SOLVER_IP_SELECT_SIMPLE 0
    int nwant;
    double pruning_tol; // >= (1.0 - SOLVER_IP_BAC_PRICE_MAXPENALTY)
} ip_exact_bac_param;

typedef struct ip_exact_bac_stats
{
    stats_item *total;
    stats_item *sep_loop;
    stats_item *sparse_edge_check;
    stats_item *full_edge_check;
    stats_item *addcuts;
    stats_item *agecuts;
    stats_item *ageedges;
    stats_item *addbad;
    stats_item *xheur;
    stats_item *misc;
    stats_item *branch_node;
    stats_item *find_branch;
    stats_item *exec_branch;
    stats_item *exec_unbranch;
    stats_item *pricing_loop;
    stats_item *recover_infeas;
    stats_item *verify_infeas_cut;
    stats_item *verify_branch_prune;
    stats_item *verify_branch_infeas;
} ip_exact_bac_stats;

typedef struct ip_branch
{
    int depth;
    int rhs;
    graph_arc *edge;
    char sense;
    int lprow;
} ip_branch;

typedef struct lp_data lp_data;
typedef struct ip_exact_bac_env
{
    int verbosity;
    ip_exact_bac_param *param;
    ip_exact_bac_stats *stats;
    lp_prob *lp;
    void *orig_prob;
    void *orig_env;
    double upperboundN;
    double lowerboundN;

    double infeas_bound;
    int (*build_lp_data)(ip_prob *ip, int col_start, int col_end,
                         lp_data **lp_data);
    int (*sep_loop)(void *prob, void *env, int *nadded);
    int (*pricing_loop)(void *prob, void *env, int *nadded);
    int (*recover_infeas)(void *prob, void *env);
    int (*verify_prune)(void *prob, void *env, int *prune);
    int (*verify_infeas)(void *prob, void *env, int *prune);
    int (*dual_bound)(void *prob, void *env, mpf_t bound);
    int (*xheur)(void *prob, void *env, void *sol, int *nadded);
    int (*check_sol)(ip_sol *sol);
    int (*get_sol)(void *prob, void *env, void *sol);
    int (*update_sol)(void *prob, void *env, int *nadded);
    int (*find_branch)(void *prob, void *env, ip_branch **branch);
    int ngot;
    int depth;
    int farkas_pricing;
    const char *branch_log;
    int branch_count;
    int history_space;
    int history_depth;
    ip_branch **history;
    int id;
    int parent_id;
    int root;
} ip_exact_bac_env;

ip_exact_bac_stats *
ip_create_exact_bac_stats(void);
void
ip_free_exact_bac_stats(ip_exact_bac_stats **stats);
int
ip_write_exact_bac_stats(ip_exact_bac_env *env);

ip_exact_bac_param *
ip_create_exact_bac_param(void);
void
ip_free_exact_bac_param(ip_exact_bac_param **param);

ip_exact_bac_env *
ip_create_exact_bac_env(void);
void
ip_free_exact_bac_env(ip_exact_bac_env **env);

ip_branch *
ip_create_branch(void);
void
ip_free_branch(ip_branch **b);
int
ip_exec_branch(ip_prob *prob, ip_exact_bac_env *env, ip_branch *b),
ip_exec_unbranch(ip_prob *prob, ip_exact_bac_env *env);

void
ip_print_branch_history(ip_exact_bac_env *env);

int
ip_branch_dfs(ip_prob *prob, ip_exact_bac_env *env, void *sol);

#endif
