#include "lp/lp.h"
#include "op-solver.h"
#include <ilcplex/cplex.h>

#undef SOLVER_CPLEX_DISPLAY

struct cplex_param
{
    int scrind;
    int simdisplay;
    int fastmip;
    int advind;
    int dpriind;
    int ppriind;
    double epper;
    double epopt;
    double eprhs;
    int perind;
    int preind;
    int aggind;
    int threads;
};

struct lp_solver
{
    CPXENVptr env;
    CPXLPptr lp;
    struct cplex_param *param;
};

lp_solver *
lp_create_solver(void)
{
    int rval = 0;
    struct lp_solver *solver;

    solver = malloc(sizeof(lp_solver));
    check_null(solver, "Out of memory", cleanup);

    solver->env = NULL;
    solver->lp  = NULL;

    solver->env = CPXopenCPLEX(&rval);
    if (rval)
    {
        fprintf(stderr, "CPXopenCPLEX failed, return code %d\n", rval);
        goto cleanup;
    }

    solver->param             = malloc(sizeof(struct cplex_param));
    solver->param->scrind     = 0;
    solver->param->simdisplay = 0;
    solver->param->fastmip    = 1;
    solver->param->advind     = 1;
    solver->param->dpriind    = CPX_DPRIIND_STEEP;
    solver->param->ppriind    = CPX_PPRIIND_STEEP;
    solver->param->epper      = 1.0E-6;
    solver->param->epopt      = 1.0E-9;
    solver->param->eprhs      = 1.0E-9;

    solver->param->perind = 0;
    solver->param->preind = 1;
    solver->param->aggind = 1;

    solver->param->threads = 1;

    return solver;

cleanup:

    return NULL;
}

void
lp_free_solver(lp_solver **solver)
{
    if (*solver)
    {
        if ((*solver)->env)
        {
            CPXfreeprob((*solver)->env, &((*solver)->lp));
            CPXcloseCPLEX(&((*solver)->env));
        }
        free((*solver)->param);
        free(*solver);
    }
}

int
lp_erase_solver(struct lp_solver *solver)
{
    int rval;

    lp_free_solver(&solver);

    solver = malloc(sizeof(lp_solver));
    check_null(solver, "Out of memory in lp_erase", cleanup);

    solver->env = NULL;
    solver->lp  = NULL;

    solver->env = CPXopenCPLEX(&rval);
    if (rval)
    {
        fprintf(stderr, "CPXopenCPLEX failed, return code %d\n", rval);
        goto cleanup;
    }

    solver->param->scrind     = 0;
    solver->param->simdisplay = 0;
    solver->param->fastmip    = 1;
    solver->param->advind     = 1;
    solver->param->dpriind    = CPX_DPRIIND_STEEP;
    solver->param->ppriind    = CPX_PPRIIND_STEEP;
    solver->param->epper      = 1.0E-6;
    solver->param->epopt      = 1.0E-9;
    solver->param->eprhs      = 1.0E-9;

    solver->param->perind = 0;
    solver->param->preind = 1;
    solver->param->aggind = 1;

cleanup:

    if (rval)
        lp_free_solver(&solver);

    return rval;
}

int
lp_change_sense(lp_prob *lp, int row, char sense)
{
    int rval = 0;
    int xindex[1];
    char asense[1];

    xindex[0] = row;
    asense[0] = sense;

    rval = CPXchgsense(lp->solver->env, lp->solver->lp, 1, xindex, asense);
    check_rval(rval, "failed", cleanup);

cleanup:
    return rval;
}

int
lp_opt(lp_prob *lp, int method)
{
    int rval = 0;

    lp->status = SOLVER_LP_UNKNOWN;
    switch (method)
    {
    case SOLVER_LP_METHOD_PRIMAL:
        rval = lp_opt_primal(lp);
        break;
    case SOLVER_LP_METHOD_DUAL:
        rval = lp_opt_dual(lp);
        break;
    default:
        rval = 1;
        fprintf(stderr, "Nonexistent method in lp_opt\n");
        break;
    }
    return rval;
}

int
lp_opt_primal(lp_prob *lp)
{
    int rval;
    int solstat;

    lp->status = SOLVER_LP_UNKNOWN;
    rval       = CPXprimopt(lp->solver->env, lp->solver->lp);
    if (rval)
    {
        if (rval == CPX_STAT_INForUNBD)
        {
            int old, oldagg;
            printf("Cplex presolve failed, switch to simplex\n");
            if (CPXgetintparam(lp->solver->env, CPX_PARAM_PREIND, &old))
            {
                fprintf(stderr, "CPXgetintparam CPX_PARAM_PREIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            if (CPXsetintparam(lp->solver->env, CPX_PARAM_PREIND, CPX_OFF))
            {
                fprintf(stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            if (CPXgetintparam(lp->solver->env, CPX_PARAM_AGGIND, &oldagg))
            {
                fprintf(stderr, "CPXgetintparam CPX_PARAM_AGGIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            if (CPXsetintparam(lp->solver->env, CPX_PARAM_AGGIND, 0))
            {
                fprintf(stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            rval = CPXprimopt(lp->solver->env, lp->solver->lp);
            if (rval)
            {
                fprintf(stderr, "CPXprimopt failed, return code %d\n", rval);
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            if (CPXsetintparam(lp->solver->env, CPX_PARAM_PREIND, old))
            {
                fprintf(stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            if (CPXsetintparam(lp->solver->env, CPX_PARAM_AGGIND, oldagg))
            {
                fprintf(stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
        }
        else
        {
            fprintf(stderr, "CPXprimopt failed, return code %d\n", rval);
            lp->status = SOLVER_LP_FAILURE;
            return SOLVER_LP_FAILURE;
        }
    }
    solstat = CPXgetstat(lp->solver->env, lp->solver->lp);
    if (solstat == CPX_STAT_INFEASIBLE)
    {
        lp->status = SOLVER_LP_INFEASIBLE;
        return SOLVER_LP_INFEASIBLE;
    }
    else if (solstat == CPX_STAT_UNBOUNDED)
    {
        lp->status = SOLVER_LP_UNBOUNDED;
        return SOLVER_LP_UNBOUNDED;
    }
    else if (solstat != CPX_STAT_OPTIMAL && solstat != CPX_STAT_OPTIMAL_INFEAS)
    {
        fprintf(stderr, "Cplex optimization status %d\n", solstat);
        lp->status = SOLVER_LP_FAILURE;
        return SOLVER_LP_FAILURE;
    }
    lp->status = SOLVER_LP_SUCCESS;
    return SOLVER_LP_SUCCESS;
}

int
lp_opt_dual(lp_prob *lp)
{
    int rval;
    int solstat;

    lp->status = SOLVER_LP_UNKNOWN;
    rval       = CPXdualopt(lp->solver->env, lp->solver->lp);

    if (rval)
    {
        if (rval == CPX_STAT_INForUNBD)
        {
            int old, oldagg;
            printf("Cplex presolve failed, switch to simplex\n");
            if (CPXgetintparam(lp->solver->env, CPX_PARAM_PREIND, &old))
            {
                fprintf(stderr, "CPXgetintparam CPX_PARAM_PREIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            if (CPXsetintparam(lp->solver->env, CPX_PARAM_PREIND, CPX_OFF))
            {
                fprintf(stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            if (CPXgetintparam(lp->solver->env, CPX_PARAM_AGGIND, &oldagg))
            {
                fprintf(stderr, "CPXgetintparam CPX_PARAM_AGGIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            if (CPXsetintparam(lp->solver->env, CPX_PARAM_AGGIND, 0))
            {
                fprintf(stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            rval = CPXdualopt(lp->solver->env, lp->solver->lp);
            if (rval)
            {
                fprintf(stderr, "CPXdualopt failed, return code %d\n", rval);
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            if (CPXsetintparam(lp->solver->env, CPX_PARAM_PREIND, old))
            {
                fprintf(stderr, "CPXsetintparam CPX_PARAM_PREIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
            if (CPXsetintparam(lp->solver->env, CPX_PARAM_AGGIND, oldagg))
            {
                fprintf(stderr, "CPXsetintparam CPX_PARAM_AGGIND failed\n");
                lp->status = SOLVER_LP_FAILURE;
                return SOLVER_LP_FAILURE;
            }
        }
        else
        {
            fprintf(stderr, "CPXdualopt failed, return code %d\n", rval);
            lp->status = SOLVER_LP_FAILURE;
            return SOLVER_LP_FAILURE;
        }
    }

    solstat = CPXgetstat(lp->solver->env, lp->solver->lp);
    if (solstat == CPX_STAT_INFEASIBLE)
    {
        lp->status = SOLVER_LP_INFEASIBLE;
        return 2;
    }
    else if (solstat != CPX_STAT_OPTIMAL && solstat != CPX_STAT_OPTIMAL_INFEAS)
    {
        fprintf(stderr, "Cplex optimization status %d\n", solstat);
        if (solstat == CPX_STAT_ABORT_IT_LIM)
        {
            int itlim;
            rval = CPXgetintparam(lp->solver->env, CPX_PARAM_ITLIM, &itlim);
            if (!rval)
            {
                printf("cplex iteration limit: %d\n", itlim);
            }
        }
        lp->status = SOLVER_LP_FAILURE;
        return SOLVER_LP_FAILURE;
    }
    lp->status = SOLVER_LP_SUCCESS;
    return 0;
}

int
lp_get_info(lp_prob *lp, lp_info **info)
{
    int rval = 0;

    *info = malloc(sizeof(lp_info));
    check_null(*info, "Out of memory in lp_get_info", cleanup);

    (*info)->ccount = 0;
    (*info)->rcount = 0;
    (*info)->cstat  = NULL;
    (*info)->rstat  = NULL;

    (*info)->ccount = CPXgetnumcols(lp->solver->env, lp->solver->lp);
    if ((*info)->ccount == 0)
    {
        fprintf(stderr, "No columns in lp_get_info\n");
        rval = 1;
        goto cleanup;
    }
    (*info)->rcount = CPXgetnumrows(lp->solver->env, lp->solver->lp);
    if ((*info)->rcount == 0)
    {
        fprintf(stderr, "No rows in lp_get_info\n");
        rval = 1;
        goto cleanup;
    }

    (*info)->cstat = malloc((*info)->ccount * sizeof(int));
    check_null((*info)->cstat, "Out of memory in lp_get_info", cleanup);
    (*info)->rstat = malloc((*info)->rcount * sizeof(int));
    check_null((*info)->rstat, "Out of memory in lp_get_info", cleanup);

    rval =
    CPXgetbase(lp->solver->env, lp->solver->lp, (*info)->cstat, (*info)->rstat);
    if (rval)
    {
        fprintf(stderr, "CPXgetbase failed\n");
        goto cleanup;
    }

    return 0;

cleanup:

    lp_free_info(info);
    return rval;
}

lp_info *
lp_create_info(int rcount, int ccount)
{
    int i;
    lp_info *info;

    info = malloc(sizeof(lp_info));
    if (info == NULL)
    {
        fprintf(stderr, "Out of memory in lp_create_info\n");
        return NULL;
    }

    info->ccount = 0;
    info->rcount = 0;
    info->cstat  = NULL;
    info->rstat  = NULL;

    if (ccount == 0)
    {
        fprintf(stderr, "No columns in lp_create_info\n");
        return NULL;
    }
    info->ccount = ccount;

    if (rcount == 0)
    {
        fprintf(stderr, "No rows in lp_create_info\n");
        return NULL;
    }
    info->rcount = rcount;

    info->cstat = malloc(info->ccount * sizeof(int));
    info->rstat = malloc(info->rcount * sizeof(int));
    if (!info->cstat || !info->rstat)
    {
        fprintf(stderr, "out of memory in lp_create_info\n");
        lp_free_info(&info);
        return NULL;
    }

    for (i = 0; i < ccount; i++)
    {
        info->cstat[i] = 0;
    }
    for (i = 0; i < rcount; i++)
    {
        info->rstat[i] = 0;
    }

    return info;
}

int
lp_is_col_active(lp_info *info, int c)
{
    if (c < 0 || c >= info->ccount)
        return 0;
    return info->cstat[c] == 1 || info->cstat[c] == 2;
}

int
lp_is_row_active(lp_info *info, int r)
{
    if (r < 0 || r >= info->rcount)
        return 0;
    return info->rstat[r] == 0;
}

void
lp_set_col_active(lp_info *info, int c)
{
    if (c >= 0 && c < info->ccount)
        info->cstat[c] = 1;
}

void
lp_set_col_inactive(lp_info *info, int c)
{
    if (c >= 0 && c < info->ccount)
        info->cstat[c] = 0;
}

void
lp_set_col_upper(lp_info *info, int c)
{
    if (c >= 0 && c < info->ccount)
        info->cstat[c] = 2;
}

void
lp_set_row_active(lp_info *info, int r)
{
    if (r >= 0 && r < info->rcount)
        info->rstat[r] = 0;
}

void
lp_set_row_inactive(lp_info *info, int r)
{
    if (r >= 0 && r < info->rcount)
        info->rstat[r] = 1;
}

void
lp_free_info(lp_info **info)
{
    if (*info != NULL)
    {
        if ((*info)->cstat)
            free((*info)->cstat);
        if ((*info)->rstat)
            free((*info)->rstat);
        free(*info);
    }
}

int
lp_get_x(lp_prob *lp, double *x)
{
    int rval = 0;
    int ncols;

    ncols = CPXgetnumcols(lp->solver->env, lp->solver->lp);
    if (ncols == 0)
    {
        fprintf(stderr, "No columns in LP\n");
        return 1;
    }
    rval = CPXgetx(lp->solver->env, lp->solver->lp, x, 0, ncols - 1);
    if (rval)
    {
        fprintf(stderr, "CPXgetx failed\n");
        return rval;
    }
    return 0;
}

int
lp_get_rc(lp_prob *lp, double *rc, int start, int end)
{
    int rval = 0;
    int ncols;

    ncols = CPXgetnumcols(lp->solver->env, lp->solver->lp);
    if (ncols == 0)
    {
        printf("WARNING: No columns in LP\n");
        return 1;
    }
    if (end + 1 > ncols)
    {
        printf("WARNING: No columns in LP\n");
        return 1;
    }

    rval = CPXgetdj(lp->solver->env, lp->solver->lp, rc, start, end);
    if (rval)
    {
        fprintf(stderr, "CPXgetdj failed\n");
        return rval;
    }
    return 0;
}

int
lp_get_slack(lp_prob *lp, double *slack)
{
    int rval = 0;
    int nrows;

    nrows = CPXgetnumrows(lp->solver->env, lp->solver->lp);
    if (nrows == 0)
    {
        fprintf(stderr, "No rows in LP\n");
        return 1;
    }
    rval = CPXgetslack(lp->solver->env, lp->solver->lp, slack, 0, nrows - 1);
    if (rval)
    {
        fprintf(stderr, "CPXgetslack failed\n");
        return rval;
    }
    return 0;
}

#if 0
static int get_farkas_multipliers(lp_prob *lp, double *y)
{
    int rval = 0;

    int i = 0, nrows, idiv, jdiv;
    double val, lb, ub;
    int *bhead = NULL;
    char *sense = NULL;

#if 1
    if (set_pareters(lp->solver->env, lp->solver->param)) {
        fprintf(stderr, "Unable to set optimization pareters\n");
    }
#endif

    if (lp->solver->env == NULL || lp->solver->lp == NULL) {
        rval = 1;
        fprintf(stderr, "env object or lp object is NULL\n");
        goto CLEANUP;
    }

    if (CPXgetmethod(lp->solver->env, lp->solver->lp) != CPX_ALG_DUAL ||
        CPXgetstat(lp->solver->env, lp->solver->lp) != CPX_STAT_INFEASIBLE) {
        rval = 1;
        fprintf(stderr, "Incorrect solution type\n");
        goto CLEANUP;
    }

    if (CPXgetijdiv(lp->solver->env, lp->solver->lp, &idiv, &jdiv)) {
        rval = 1;
        fprintf(stderr, "CPXgetijdiv failed\n");
        goto CLEANUP;
    }

    if ((jdiv == -1 && idiv == -1) || (jdiv != -1 && idiv != -1)) {
        rval = 1;
        fprintf(stderr, "CPLEX returned illegal indices\n");
        goto CLEANUP;
    }

    nrows = CPXgetnumrows(lp->solver->env, lp->solver->lp);
    if (nrows == 0) {
        rval = 1;
        fprintf(stderr, "solver->lp has no rows\n");
        goto CLEANUP;
    }

    bhead = malloc(nrows * sizeof(int));
    sense = malloc(nrows * sizeof(char));
    if (bhead == NULL || sense == NULL) {
        rval = -1;
        fprintf(stderr, "Out of memory\n");
        goto CLEANUP;
    }

    if (CPXgetbhead(lp->solver->env, lp->solver->lp, bhead, NULL)) {
        rval = 1;
        fprintf(stderr, "CPXgetbhead failed\n");
        goto CLEANUP;
    }

    if (CPXgetsense(lp->solver->env, lp->solver->lp, sense, 0, nrows - 1)) {
        rval = 1;
        fprintf(stderr, "CPXgetsense failed\n");
        goto CLEANUP;
    }

    if (jdiv >= 0) {
        for (i = 0; i < nrows; i++) {
            if (bhead[i] == jdiv)
                break;
        }
        if (i == nrows) {
            rval = 1;
            fprintf(stderr, "Basis index not found\n");
            goto CLEANUP;
        }
        if (CPXgetx(lp->solver->env, lp->solver->lp, &val, jdiv, jdiv)) {
            rval = 1;
            fprintf(stderr, "CPXgetx failed\n");
            goto CLEANUP;
        }
        if (CPXgetlb(lp->solver->env, lp->solver->lp, &lb, jdiv, jdiv)) {
            rval = 1;
            fprintf(stderr, "CPXgetlb failed\n");
            goto CLEANUP;
        }
        if (CPXgetub(lp->solver->env, lp->solver->lp, &ub, jdiv, jdiv)) {
            rval = 1;
            fprintf(stderr, "CPXgetub failed\n");
            goto CLEANUP;
        }
    } else {
        for (i = 0; i < nrows; i++) {
            if (bhead[i] == -idiv - 1)
                break;
        }
        if (i == nrows) {
            rval = 1;
            fprintf(stderr, "Basis index not found\n");
            goto CLEANUP;
        }
        if (CPXgetslack(lp->solver->env, lp->solver->lp, &val, idiv, idiv)) {
            rval = 1;
            fprintf(stderr, "CPXgetslack failed\n");
            goto CLEANUP;
        }
        lb = 0.0;
        if (sense[idiv] == 'E')
            ub = 0.0;
        else
            ub = CPX_INFBOUND;
        if (sense[idiv] == 'G')
            val *= -1.0;
    }

    if (CPXbinvrow(lp->solver->env, lp->solver->lp, i, y)) {
        rval = 1;
        fprintf(stderr, "CPXbinvrow failed\n");
        goto CLEANUP;
    }

    if (val < lb)
        for (i = 0; i < nrows; i++)
            y[i] *= -1.0;

    for (i = 0; i < nrows; i++) {
        if (sense[i] == 'L' && y[i] > 0.0)
            y[i] = 0.0;
        if (sense[i] == 'G' && y[i] < 0.0)
            y[i] = 0.0;
    }

CLEANUP:

    if (bhead)
        free(bhead);
    if (sense)
        free(sense);

    return rval;
}
#endif

int
lp_get_pi(lp_prob *lp, double *pi)
{
    int rval = 0;
    int nrows;

#if 0
  if ( CPXgetmethod (lp->solver->env, lp->solver->lp) == CPX_ALG_DUAL       &&
       CPXgetstat (lp->solver->env, lp->solver->lp )  == CPX_STAT_INFEASIBLE  )
  { rval = get_farkas_multipliers (lp->solver, pi);
    printf(" >FARKAS\n");
    if (rval)
    { fprintf (stderr, "get_farkas_multipliers failed\n"); return rval;
    }
    return 2;
  }
#else
    if (CPXgetmethod(lp->solver->env, lp->solver->lp) == CPX_ALG_DUAL &&
        CPXgetstat(lp->solver->env, lp->solver->lp) == CPX_STAT_INFEASIBLE)
    {
        double proof;
        rval = CPXdualfarkas(lp->solver->env, lp->solver->lp, pi, &proof);
        if (rval)
        {
            fprintf(stderr, "get_farkas_multipliers failed\n");
            return rval;
        }
        return 2;
    }
#endif

    nrows = CPXgetnumrows(lp->solver->env, lp->solver->lp);
    if (nrows == 0)
    {
        fprintf(stderr, "No rows in LP\n");
        return 1;
    }
    rval = CPXgetpi(lp->solver->env, lp->solver->lp, pi, 0, nrows - 1);
    if (rval)
    {
        fprintf(stderr, "CPXgetpi failed\n");
        return rval;
    }
    return 0;
}

double
lp_get_objval(lp_prob *lp)
{
    int rval;
    double obj;

    rval = CPXgetobjval(lp->solver->env, lp->solver->lp, &obj);
    if (rval)
    {
        fprintf(stderr, "CPXgetobjval failed\n");
        return -rval;
    }
    return obj;
}

int
lp_get_nrows(lp_prob *lp)
{
    return CPXgetnumrows(lp->solver->env, lp->solver->lp);
}

int
lp_get_ncols(lp_prob *lp)
{
    return CPXgetnumcols(lp->solver->env, lp->solver->lp);
}

int
lp_get_nnonzeros(lp_prob *lp)
{
    return CPXgetnumnz(lp->solver->env, lp->solver->lp);
}

int
lp_add_rows(lp_prob *lp, lp_data *data)
{
    int rval = 0;

    rval = CPXaddrows(lp->solver->env, lp->solver->lp, 0, data->nrows,
                      data->nzcnt, data->rhs, data->sense, data->beg, data->ind,
                      data->val, NULL, NULL);
    if (rval)
        fprintf(stderr, "CPXaddrows failed\n");
    return rval;
}

int
lp_add_cols(lp_prob *lp, lp_data *data)
{
    int rval = 0;

    rval = CPXaddcols(lp->solver->env, lp->solver->lp, data->ncols, data->nzcnt,
                      data->obj, data->beg, data->ind, data->val, data->lb,
                      data->ub, NULL);
    if (rval)
        fprintf(stderr, "CPXaddcols failed\n");
    return rval;
}

int
lp_del_row(lp_prob *lp, int i)
{
    int rval = 0;
    int locali[1];

    locali[0] = i;
    if (CPXpivotin(lp->solver->env, lp->solver->lp, locali, 1))
    {
        fprintf(stderr, "CPXpivotin failed, continuing anyway\n");
    }
    rval = CPXdelrows(lp->solver->env, lp->solver->lp, i, i);
    if (rval)
        fprintf(stderr, "CPXdelrows failed\n");
    return rval;
}

int
lp_del_rows(lp_prob *lp, int *delstat)
{
    int rval     = 0;
    int *dellist = NULL;
    int delcnt   = 0;
    int i;
    int j;
    int rcnt = CPXgetnumrows(lp->solver->env, lp->solver->lp);

    for (i = 0; i < rcnt; i++)
    {
        if (delstat[i])
            delcnt++;
    }
    if (delcnt == 0)
    {
        printf("WARNING: lp_del_rows with no deleted rows\n");
        return 0;
    }

    dellist = malloc(delcnt * sizeof(int));
    check_null(dellist, "Out of memory in lp_del_rows", cleanup);

    for (i = 0, j = 0; i < rcnt; i++)
    {
        if (delstat[i])
            dellist[j++] = i;
    }
    if (j != delcnt)
    {
        printf("Lost some deleted rows\n");
        free(dellist);
        return 1;
    }

    if (CPXpivotin(lp->solver->env, lp->solver->lp, dellist, delcnt))
        fprintf(stderr, "WARNING: CPXpivotin failed, continuing anyway\n");

    rval = CPXdelsetrows(lp->solver->env, lp->solver->lp, delstat);
    check_rval(rval, "CPXdelsetrows failed", cleanup);

cleanup:
    if (dellist)
        free(dellist);

    return rval;
}

int
lp_del_col(lp_prob *lp, int i)
{
    int rval = 0;
    int locali[1];
    char lu[1];
    double bd[1];

    locali[0] = i;
    lu[0]     = 'B';
    bd[0]     = 0.0;

    if (CPXchgbds(lp->solver->env, lp->solver->lp, 1, locali, lu, bd))
    {
        printf("WARNING: CPXchgbds failed, continuing anyway\n");
    }

    if (CPXdualopt(lp->solver->env, lp->solver->lp))
    {
        printf("WARNING: CPXdualopt failed, continuing anyway\n");
    }

    if (CPXpivotout(lp->solver->env, lp->solver->lp, locali, 1))
    {
        printf("WARNING: CPXpivotout failed, continuing anyway\n");
    }

    rval = CPXdelcols(lp->solver->env, lp->solver->lp, i, i);
    if (rval)
        fprintf(stderr, "CPXdelcols failed\n");
    return rval;
}

int
lp_del_cols(lp_prob *lp, int *delstat)
{
    int rval     = 0;
    int *dellist = NULL;
    char *lu     = NULL;
    double *bd   = NULL;
    int delcnt   = 0;
    int i;
    int j;
    int ccnt = CPXgetnumcols(lp->solver->env, lp->solver->lp);

    for (i = 0; i < ccnt; i++)
    {
        if (delstat[i])
            delcnt++;
    }
    if (delcnt == 0)
    {
        printf("WARNING: delete_set_of_columns with no deleted columns\n");
        return 0;
    }
    dellist = malloc(delcnt * sizeof(int));
    lu      = malloc(delcnt * sizeof(char));
    bd      = malloc(delcnt * sizeof(double));
    if (dellist == NULL || lu == NULL || bd == NULL)
    {
        printf("WARNING: Out of memory in delete_set_of_columns\n");
        free(dellist);
        free(lu);
        free(bd);
        return 1;
    }
    for (i = 0, j = 0; i < ccnt; i++)
    {
        if (delstat[i])
        {
            lu[j]        = 'B';
            bd[j]        = 0.0;
            dellist[j++] = i;
        }
    }
    if (j != delcnt)
    {
        printf("WARNING: Lost some deleted columns\n");
        free(dellist);
        free(lu);
        free(bd);
        return 1;
    }

    if (CPXchgbds(lp->solver->env, lp->solver->lp, delcnt, dellist, lu, bd))
        printf("WARNING: CPXchgbds failed, stumbling on anyway\n");

    if (CPXdualopt(lp->solver->env, lp->solver->lp))
        printf("WARNING: CPXdualopt failed, continuing anyway\n");

    if (CPXpivotout(lp->solver->env, lp->solver->lp, dellist, delcnt))
        printf("WARNING: CPXpivotout failed, continuing anyway\n");

    free(dellist);
    free(lu);
    free(bd);

    rval = CPXdelsetcols(lp->solver->env, lp->solver->lp, delstat);
    if (rval)
        printf("WARNING: CPXdelsetcols failed\n");
    return rval;
}

int
lp_set_bnd(lp_prob *lp, int col, char lower_or_upper, double bnd)
{
    int cindex[1];
    double bd[1];
    char lu[1];
    int rval;

    cindex[0] = col;
    lu[0]     = lower_or_upper;
    bd[0]     = bnd;

    rval = CPXchgbds(lp->solver->env, lp->solver->lp, 1, cindex, lu, bd);
    if (rval)
    {
        fprintf(stderr, "Couldn't set bnd on variable %d in cplex\n", col);
        return rval;
    }
    return 0;
}

int
lp_init(lp_prob *lp, lp_data *data)
{
    int rval = 0;

    lp->solver->lp = CPXcreateprob(lp->solver->env, &rval, "probname");
    if (!lp->solver->lp || rval)
    {
        fprintf(stderr, "CPXcreateprob failed, return code %d\n", rval);
        return 1;
    }

    rval =
    CPXcopylp(lp->solver->env, lp->solver->lp, data->ncols, data->nrows,
              data->objsense, data->obj, data->rhs, data->sense, data->beg,
              data->cnt, data->ind, data->val, data->lb, data->ub, NULL);
    if (rval)
    {
        fprintf(stderr, "CPXcopylp failed, return code %d\n", rval);
        return 1;
    }

    return 0;
}
