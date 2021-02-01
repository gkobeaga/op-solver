#include "kp/kp.h"
#include "op-solver.h"

static void
kp_init_prob_work(kp_prob *kp, int n)
{
    int i;
    kp->n = n;
    kp->p = malloc(n * sizeof(double));
    memset(kp->p, 0.0, n * sizeof(double));
    kp->w = malloc(n * sizeof(double));
    memset(kp->w, 0.0, n * sizeof(double));
    kp->tot_profit = 0.0;
    kp->cap        = 0.0;
    kp->ord_ind    = malloc(n * sizeof(int));
    for (i = 0; i < n; i++) kp->ord_ind[i] = i;
    kp->ord_w    = NULL;
    kp->ord_p    = NULL;
    kp->min_w    = NULL;
    kp->sol_stat = SOLVER_UNDEF;
    kp->sol      = kp_create_sol(kp);
    return;
}

kp_prob *
kp_create_prob(solver_data *data)
{
    kp_prob *kp = malloc(sizeof(kp_prob));
    kp_init_prob_work(kp, data->n);
    return kp;
}

static void
kp_delete_prob_work(kp_prob *kp);

void
kp_erase_prob(kp_prob *kp)
{
    int n = kp->n;
    kp_delete_prob_work(kp);
    kp_init_prob_work(kp, n);
    return;
}

static void
kp_delete_prob_work(kp_prob *kp)
{
    if (kp->w)
        free(kp->w);
    if (kp->p)
        free(kp->p);
    if (kp->min_w)
        free(kp->min_w);
    if (kp->ord_ind)
        free(kp->ord_ind);
    if (kp->ord_p)
        free(kp->ord_p);
    if (kp->ord_w)
        free(kp->ord_w);
    kp_free_sol(&(kp->sol));
    return;
}

void
kp_free_prob(kp_prob **kp)
{
    if (*kp)
    {
        kp_delete_prob_work(*kp);
        free(*kp);
        *kp = NULL;
    }
    return;
}
