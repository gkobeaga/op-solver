#include "kp/kp.h"
#include "op-solver.h"

static void
kp_create_sol_work(kp_sol *sol, int n)
{
    sol->tot_n     = n;
    sol->selected  = malloc(n * sizeof(int));
    sol->sposition = malloc(n * sizeof(int));
    sol->val       = 0.0;
    sol->weight    = 0.0;
    sol->ns        = 0;
    memset(sol->selected, 0, n * sizeof(int));
    memset(sol->sposition, n, n * sizeof(int));
}

kp_sol *
kp_create_sol(kp_prob *kp)
{
    kp_sol *sol = malloc(sizeof(kp_sol));
    kp_create_sol_work(sol, kp->n);

    return sol;
}

static void
kp_delete_sol_work(kp_sol *sol)
{
    free(sol->selected);
    free(sol->sposition);
}

void
kp_erase_sol(kp_sol *sol)
{
    int n = sol->tot_n;
    kp_delete_sol_work(sol);
    kp_create_sol_work(sol, n);
    return;
}

void
kp_free_sol(kp_sol **sol)
{
    kp_delete_sol_work(*sol);
    free(*sol);
    *sol = NULL;
    return;
}

int
kp_start_sol(kp_prob *kp, kp_env *env, kp_sol *sol)
{
    int rval = 0;
    int i;
    double residual_cap, tmpbound, tmpprofit;
    kp_param *param = env->param;

    check_null(sol, "kp_start_solution failed: sol not created", CLEANUP);

    kp_erase_sol(sol);
    residual_cap = kp->cap;
    for (i = 0; i < kp->n && kp->w[i] <= residual_cap; i++)
    {
        sol->val += kp->p[i];
        residual_cap -= kp->w[i];
        sol->weight += kp->w[i];
        sol->selected[kp->ord_ind[i]] = 1;
        sol->ns++;
    }

    i--;
    if (i < kp->n - 2)
    {
        if (env->p_integer)
            tmpprofit =
            ceil((kp->w[i + 1] - residual_cap) * kp->p[i] / kp->w[i]);
        else
            tmpprofit = (kp->w[i + 1] - residual_cap) * kp->p[i] / kp->w[i] +
                        sol->val * param->epsilon;
        kp->upperbound = sol->val + kp->p[i + 1] - tmpprofit;

        if (env->p_integer)
            tmpprofit = floor(residual_cap * kp->p[i + 2] / kp->w[i + 2]);
        else
            tmpprofit = residual_cap * kp->p[i + 2] / kp->w[i + 2] -
                        sol->val * param->epsilon;

        tmpbound = sol->val + tmpprofit;
        if (kp->upperbound < tmpbound)
            kp->upperbound = tmpbound;
    }
    else if (i == kp->n - 2)
    {
        if (env->p_integer)
            tmpprofit =
            ceil((kp->w[i + 1] - residual_cap) * kp->p[i] / kp->w[i]);
        else
            tmpprofit = (kp->w[i + 1] - residual_cap) * kp->p[i] / kp->w[i] +
                        sol->val * param->epsilon;
        kp->upperbound = sol->val + kp->p[i + 1] - tmpprofit;
    }
    else if (i == kp->n - 1)
        kp->upperbound = sol->val;

CLEANUP:
    return rval;
}
