#ifndef LP_H
#define LP_H

#include "../../op-solver.h"

#define SOLVER_LP_METHOD_DUAL 1
#define SOLVER_LP_METHOD_PRIMAL 2
#define SOLVER_LP_METHOD_BARRIER 3

#define SOLVER_LP_SUCCESS 0
#define SOLVER_LP_FAILURE 1
#define SOLVER_LP_UNBOUNDED 2
#define SOLVER_LP_INFEASIBLE 3
#define SOLVER_LP_UNKNOWN 4

/* used to check if lp solis integral  */
#define SOLVER_INTTOL (0.000001)

typedef struct lp_sol lp_sol;
typedef struct lp_solver lp_solver;
typedef struct lp_info lp_info;
typedef void lp_cut;

typedef struct lp_prob
{
    solver_graph *graph;
    lp_solver *solver;
    int sense;
    int infeasible;
    int status;
    lp_sol *sol;

} lp_prob;

typedef struct lp_param
{
    int init_phase;
} lp_param;

typedef struct lp_stats
{
    stats_item *total;
} lp_stats;

typedef struct lp_env
{
    int verbosity;
    lp_param *param;
    lp_stats *stats;
} lp_env;

typedef struct lp_data
{
    const char *name;
    int objsense;
    double *obj;
    double *rhs;
    char *sense;
    int *beg;
    int *cnt;
    int *ind;
    double *val;
    double *lb;
    double *ub;
    int nrows;
    int ncols;
    int nzcnt;
    int rowspace;
    int colspace;
    int nzspace;
} lp_data;

typedef struct lp_sol
{
    double val;
    solver_graph *graph;
    double *x;
    double *rc;
    int integral;
} lp_sol;

typedef struct lp_info
{
    int rcount;
    int ccount;
    int *rstat;
    int *cstat;
} lp_info;

lp_prob *
lp_create_prob(void);
void
lp_free_prob(lp_prob **prob);

lp_solver *
lp_create_solver(void);
void
lp_free_solver(lp_solver **solver);
int
lp_erase_solver(lp_solver *solver);
int
lp_change_sense(lp_prob *lp, int row, char sense);
int
lp_opt(lp_prob *lp, int method);
int
lp_opt_primal(lp_prob *lp),
lp_opt_dual(lp_prob *lp);

lp_info *
lp_create_info(int rcount, int ccount);
void
lp_free_info(lp_info **info);
int
lp_get_info(lp_prob *lp, lp_info **info),
lp_is_col_active(lp_info *info, int c), lp_is_row_active(lp_info *info, int r);
void
lp_set_col_active(lp_info *info, int c),
lp_set_col_inactive(lp_info *info, int c),
lp_set_col_upper(lp_info *info, int c), lp_set_row_active(lp_info *info, int r),
lp_set_row_active(lp_info *info, int r), lp_free_info(lp_info **info);
int
lp_get_x(lp_prob *lp, double *x),
lp_get_rc(lp_prob *lp, double *rc, int start, int end),
lp_get_slack(lp_prob *lp, double *slack), lp_get_pi(lp_prob *lp, double *pi);
double
lp_get_objval(lp_prob *lp);
int
lp_get_nrows(lp_prob *lp),
lp_get_ncols(lp_prob *lp), lp_get_nnonzeros(lp_prob *lp),
lp_add_rows(lp_prob *lp, lp_data *data),
lp_add_cols(lp_prob *lp, lp_data *data), lp_del_row(lp_prob *lp, int i),
lp_del_rows(lp_prob *lp, int *delstat), lp_del_col(lp_prob *lp, int i),
lp_del_cols(lp_prob *lp, int *delstat),
lp_set_bnd(lp_prob *lp, int col, char lower_or_upper, double bnd),
lp_init(lp_prob *lp, lp_data *data),
lp_opt_limited_dualopt(lp_prob *lp, int iterationlim, int *status,
                       double *objupperlim);

lp_env *
lp_create_env(void);
void
lp_free_env(lp_env **env);
lp_param *
lp_create_param(void);
void
lp_free_param(lp_param **param);
lp_stats *
lp_create_stats(void);
void
lp_free_stats(lp_stats **stats);

lp_sol *
lp_create_sol(void);
void
lp_free_sol(lp_sol **sol);
int
lp_update_sol(lp_prob *lp, lp_sol *sol);

lp_data *
lp_create_data(int ncols, int nrows, int nzcnt);
void
lp_free_data(lp_data **data);
int
lp_realloc_data(lp_data *data, int count, double scale);

#endif
