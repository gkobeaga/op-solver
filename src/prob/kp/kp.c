#include "kp/kp.h"
#include "op-solver.h"

static int
check_input(kp_prob *kp, kp_env *env);
static int
reorder_items_begining(kp_prob *kp, kp_env *env);
static int
reorder_items_ending(kp_prob *kp, kp_env *env);

int
kp_opt(kp_prob *kp, kp_env *env, kp_sol *sol)
{
    int rval        = 0;
    kp_stats *stats = env->stats;
    kp_param *param = env->param;

    rval = stats_start(stats->total);
    check_rval(rval, "stats_start failed", cleanup);

    /**************************************************************************/
    /* Check input data                                                       */
    /**************************************************************************/
    if (param->check_input == 1)
    {
        rval = check_input(kp, env);
        if (rval == 5)
        {
            rval                 = 0;
            param->reorder_items = 1;
        }
        else if (rval != 0)
            printf("kp rval %d\n", rval);
        check_rval(rval, "Wrong input data in the Knapsack Problem\n", done);
    }

    if (param->reorder_items == 1)
        reorder_items_begining(kp, env);

    /**************************************************************************/
    /* Initial solution                                                       */
    /**************************************************************************/
    if (env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("\n");
        printf("kp   : Building initial solution.\n");
    }

    rval = stats_start(stats->init_sol);
    check_rval(rval, "stats_start failed", cleanup);
    if (!(sol))
        sol = kp_create_sol(kp);
    kp_start_sol(kp, env, sol);
    rval = stats_stop(stats->init_sol, 1);
    check_rval(rval, "stats_stop failed", cleanup);

    if (env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("kp   : Value: %.2f , Weight: %.2f, Visited: %d\n", sol->val,
               sol->weight, sol->ns);
        printf("kp   : Upper bound: %f\n", kp->upperbound);
        printf("kp   : Time: %.2ld sec \n",
               stats_get_current_time(stats->init_sol));
    }

    /**************************************************************************/
    /* Solve exactly                                                          */
    /**************************************************************************/
    if (param->exact_tech == KP_EXACT_NONE)
    {
        if (env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("kp   : Skipping to solve exactly.\n");
    }
    else if (param->exact_tech == KP_EXACT_BAB)
    {
        if (env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("kp   : Starting branch-and-bound method...\n");
        kp_opt_bab(kp, env, sol);
    }
    else
    {
        if (env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("kp   : exact method not implemented.\n");
    }

    rval = stats_stop(stats->total, 1);
    check_rval(rval, "stats_stop failed", cleanup);

    /**************************************************************************/
    /* Print Info                                                             */
    /**************************************************************************/
    if (env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        if (kp->sol_stat == SOLVER_OPT || kp->sol_stat == SOLVER_FEAS)
        {
            printf("\n");
            printf("kp  : Best solution value: %.2f\n", sol->val);
            printf("kp  : Best solution weight: %.2f\n", sol->weight);
            printf("kp   : Time: %.2ld sec \n",
                   stats_get_current_time(stats->total));
        }
        else if (kp->sol_stat == SOLVER_NOFEAS)
        {
            printf("\n");
            printf("kp   : Error! No feasible solution found\n");
            rval = SOLVER_NOFEAS;
        }
        else
        {
            printf("\n");
            printf("kp   : Error! Problem not solved\n");
            rval = SOLVER_UNDEF;
        }
    }

    rval = EXIT_SUCCESS;
done:

    if (param->reorder_items == 1)
        reorder_items_ending(kp, env);

    return rval;

cleanup:

    rval = stats_stop(stats->total, 0);
    check_rval(rval, "stats_stop failed", cleanup);

    return rval;
}

typedef struct kp_item
{
    int ind;
    double p;
    double w;
} kp_item;

static int
sort_items(const void *xx, const void *yy)
{
    double val1, val2;
    kp_item *x = (kp_item *)xx, *y = (kp_item *)yy;

    val1 = (x->p) / x->w;
    val2 = (y->p) / y->w;

    if (val1 > val2)
        return -1;
    if (val1 < val2)
        return 1;

    return 0;
}

static int
reorder_items_begining(kp_prob *kp, kp_env *env)
{
    int rval = 0;
    int i;
    kp_item *items = NULL;

    kp->ord_w = malloc(kp->n * sizeof(double));
    check_null(kp->ord_w, "out of memory", cleanup);
    kp->ord_p = malloc(kp->n * sizeof(double));
    check_null(kp->ord_p, "out of memory", cleanup);

    items = malloc(kp->n * sizeof(kp_item));
    check_null(items, "out of memory", cleanup);

    for (i = 0; i < kp->n; i++)
    {
        items[i].ind = i;
        kp->ord_p[i] = kp->p[i];
        items[i].p   = kp->p[i];
        kp->ord_w[i] = kp->w[i];
        items[i].w   = kp->w[i];
    }

    qsort(items, kp->n, sizeof(kp_item), sort_items);

    for (i = 0; i < kp->n; i++)
    {
        kp->ord_ind[i] = items[i].ind;
        kp->p[i]       = items[i].p;
        kp->w[i]       = items[i].w;
    }

cleanup:
    if (items)
        free(items);
    return rval;
}

static int
reorder_items_ending(kp_prob *kp, kp_env *env)
{
    double p, w;
    int *ind = NULL;
    int i;

    ind = malloc(kp->n * sizeof(int));

    for (i = 0; i < kp->n; i++)
    {
        ind[kp->ord_ind[i]] = i;
        p                   = kp->p[i];
        w                   = kp->w[i];
        kp->p[i]            = kp->ord_p[i];
        kp->w[i]            = kp->ord_w[i];
        kp->ord_p[i]        = p;
        kp->ord_w[i]        = w;
    }

    free(kp->ord_ind);
    kp->ord_ind = ind;

    return 0;
}

static int
check_input(kp_prob *kp, kp_env *env)
{
    int rval = 0;
    double oldrel, rel, tot_weight;
    int i;

    if (kp->n < 2)
    {
        rval = 1;
        goto CLEANUP;
    }

    if (kp->cap <= 0.0)
    {
        rval = 2;
        goto CLEANUP;
    }

    tot_weight = 0.0;
    rel        = kp->p[0] / kp->w[0];

    for (i = 0; i < kp->n; i++)
    {
        oldrel = rel;

        if (kp->p[i] <= 0.0)
        {
            rval = 2;
            goto CLEANUP;
        }

        if (kp->w[i] <= 0.0)
        {
            rval = 2;
            goto CLEANUP;
        }

        tot_weight += kp->w[i];
        if (kp->w[i] > kp->cap)
        {
            rval = 3;
            return 0;
        }

        rel = kp->p[i] / kp->w[i];
        if (rel > oldrel)
        {
            rval = 5;
            goto CLEANUP;
        }
    }

    if (tot_weight <= kp->cap)
        rval = 4;

CLEANUP:
    return rval;
}

int
kp_get_mincover(kp_prob *kp, kp_sol *sol)
{
    int rval = 0;
    int i, ind;

    if (kp->sol_stat == SOLVER_UNDEF)
    {
        rval = 1;
        goto CLEANUP;
    }

    i = 0;
    do
    {
        ind = kp->ord_ind[i];
        if (sol->selected[ind] == 0 && sol->weight + kp->w[ind] <= kp->cap)
        {
            sol->selected[ind] = 1;
            sol->weight += kp->w[ind];
            sol->ns += 1;
        }
        i++;
    } while (i < kp->n && sol->weight + kp->min_w[ind] <= kp->cap);

CLEANUP:
    return rval;
}
