#include "lp/lp.h"
#include "op-solver.h"

struct dummy_param
{
};

struct lp_solver
{
    void *env;
    void *lp;
    struct dummy_param *param;
};

lp_solver *
lp_create_solver(void)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return NULL;
}

void
lp_free_solver(lp_solver **solver)
{
    fprintf(stderr, "Using dummy LP solver.\n");
}

int
lp_erase_solver(struct lp_solver *solver)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_change_sense(lp_prob *lp, int row, char sense)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_opt(lp_prob *lp, int method)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_opt_primal(lp_prob *lp)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_opt_dual(lp_prob *lp)
{
    fprintf(stderr, "Using dummy LP solver.\n");

    return 1;
}

int
lp_get_info(lp_prob *lp, lp_info **info)
{
    fprintf(stderr, "Using dummy LP solver.\n");

    return 1;
}

lp_info *
lp_create_info(int rcount, int ccount)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return NULL;
}

int
lp_is_col_active(lp_info *info, int c)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_is_row_active(lp_info *info, int r)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

void
lp_set_col_active(lp_info *info, int c)
{
    fprintf(stderr, "Using dummy LP solver.\n");
}

void
lp_set_col_inactive(lp_info *info, int c)
{
    fprintf(stderr, "Using dummy LP solver.\n");
}

void
lp_set_col_upper(lp_info *info, int c)
{
    fprintf(stderr, "Using dummy LP solver.\n");
}

void
lp_set_row_active(lp_info *info, int r)
{
    fprintf(stderr, "Using dummy LP solver.\n");
}

void
lp_set_row_inactive(lp_info *info, int r)
{
    fprintf(stderr, "Using dummy LP solver.\n");
}

void
lp_free_info(lp_info **info)
{
    fprintf(stderr, "Using dummy LP solver.\n");
}

int
lp_get_x(lp_prob *lp, double *x)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_get_rc(lp_prob *lp, double *rc, int start, int end)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_get_slack(lp_prob *lp, double *slack)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_get_pi(lp_prob *lp, double *pi)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

double
lp_get_objval(lp_prob *lp)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_get_nrows(lp_prob *lp)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_get_ncols(lp_prob *lp)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_get_nnonzeros(lp_prob *lp)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_add_rows(lp_prob *lp, lp_data *data)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_add_cols(lp_prob *lp, lp_data *data)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_del_row(lp_prob *lp, int i)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_del_rows(lp_prob *lp, int *delstat)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_del_col(lp_prob *lp, int i)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_del_cols(lp_prob *lp, int *delstat)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_set_bnd(lp_prob *lp, int col, char lower_or_upper, double bnd)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}

int
lp_init(lp_prob *lp, lp_data *data)
{
    fprintf(stderr, "Using dummy LP solver.\n");
    return 1;
}
