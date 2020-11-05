#ifndef CP_H
#define CP_H

#include "../../op-solver.h"

typedef struct cp_sol
{
    int tot_n;
    int *cod_fr; /* Forward codification */
    int *cod_bk; /* Backward codification */
    int *cycle;
    int *selected;
    int *sposition;
    double val;
    double cap;
    int ns;
} cp_sol;

typedef struct cp_pop
{
    int size;
    int stop_per;
    cp_sol **sol;
    int *rankperm;
    double mean_val;
    double best_val;
    int best_ind;
    double q25_val;
    int q25_ind;
    double q50_val;
    int q50_ind;
    double q75_val;
    int q75_ind;
    double stop_val;
    int stop_ind;
    double worst_val;
    int worst_ind;
    int *parent;
} cp_pop;

typedef struct ip_prob ip_prob;
typedef struct cp_prob cp_prob;
typedef struct cp_env cp_env;
typedef struct cp_prob
{
    int n;
    double cap; /* Total distance limit */
    int from;   /* Departure node */
    int to;     /* Arraiving node */
    cp_pop *pop;
    solver_data *data;
    cp_sol *sol;
    int sol_status;
    // Init
    ip_prob *ip;
    int (*eval_sol_obj)(cp_prob *cp, cp_sol *sol);
} cp_prob;

typedef struct cp_param
{
    long time_limit;
    int check_input;
    int reorder_items;
    int appr;
#define SOLVER_CP_APPR_UNDEFINED 0
#define SOLVER_CP_APPR_EXACT_BAC 1
#define SOLVER_CP_APPR_HEUR_EA 2
} cp_param;

typedef struct cp_stats
{
    stats_item *total;
} cp_stats;

typedef struct ip_env ip_env;
typedef struct cp_init_env cp_init_env;
typedef struct cp_heur_env cp_heur_env;
typedef struct cp_exact_env cp_exact_env;
typedef struct cp_env
{
    int verbosity;
    cp_param *param;
    cp_stats *stats;
    cp_init_env *init;
    cp_heur_env *heur;
    cp_exact_env *exact;
    ip_env *ip;
    char sol_file[50];
} cp_env;

cp_prob *
cp_create_prob(solver_data *data);
void
cp_erase_prob(cp_prob *cp);
void
cp_free_prob(cp_prob **cp);
int
cp_read_prob(cp_prob *cp, cp_env *env, const char *fname, int flags);

int
cp_opt(cp_prob *cp, cp_env *env, cp_sol *sol);

cp_env *
cp_create_env(void);
void
cp_free_env(cp_env **env);
int
cp_parse_args(int argc, char *argv[], cp_env *env);
cp_param *
cp_create_param(void);
void
cp_free_param(cp_param **param);
cp_stats *
cp_create_stats(void);
void
cp_free_stats(cp_stats **stats);

cp_sol *
cp_create_sol(cp_prob *cp);
void
cp_copy_sol(cp_sol *insol, cp_sol *outsol);
void
cp_erase_sol(cp_sol *sol);
void
cp_free_sol(cp_sol **sol);
int
cp_get_sol_from_graph(cp_prob *cp, solver_graph *graph, cp_sol **sol);
void
cp_print_sol(cp_prob *cp, cp_sol *sol),
cp_plot_sol(cp_prob *cp, cp_sol *sol);
int
cp_write_sol(cp_prob *cp, cp_sol *sol, const char *fname);

solver_graph *
cp_conv_sol_to_graph(cp_prob *cp, cp_sol *sol);

cp_pop *
cp_create_pop(cp_prob *cp, int size);
void
cp_set_pop_sol(cp_pop *pop, cp_sol *sol, int pos);
int
cp_update_pop(cp_pop *pop);
void
cp_delete_pop(cp_pop *pop);
void
cp_erase_pop(cp_pop *pop);
void
cp_free_pop(cp_pop **pop);

// LIB: TSP
typedef struct cp_sol tsp_sol;
void
cp_conv_sol_cp2tsp(solver_data *data, cp_sol *cpsol, tsp_sol *tspsol),
cp_conv_sol_tsp2cp(solver_data *data, tsp_sol *tspsol, cp_sol *cpsol);

#include "exact/exact.h"
#include "heur/heur.h"
#include "init/init.h"

#endif
