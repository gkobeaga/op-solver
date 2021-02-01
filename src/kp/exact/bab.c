#include "kp/kp.h"
#include "op-solver.h"

#define BIGDOUBLE (1e30)

typedef struct branch_node
{
    int level;
    int profit;
    double cap;
    double upperbound;
    int branch;
    int unbranch;
    int *x;
    int bbcount;
} branch_node;

static int
initialize(kp_prob *kp, kp_env *env, kp_sol *sol, branch_node *node);
static int
dfs(kp_prob *kp, kp_env *env, kp_sol *sol, branch_node *node);
static int
branch(kp_prob *kp, kp_env *env, kp_sol *sol, branch_node *node);
static int
unbranch(kp_prob *kp, kp_env *env, branch_node *node);

static int
save_opt_sol(kp_prob *kp, kp_env *env, kp_sol *sol, branch_node *node);

static double
get_upperbound(kp_prob *kp, kp_env *env, branch_node *node);

int
kp_opt_bab(kp_prob *kp, kp_env *env, kp_sol *sol)
{
    int rval          = 0;
    branch_node *node = NULL;

    node    = malloc(sizeof(branch_node));
    node->x = NULL;

    node->x = malloc(kp->n * sizeof(int));
    check_null(node->x, "Out of memory in kp", CLEANUP);

    kp->min_w = malloc(kp->n * sizeof(double));
    check_null(kp->min_w, "Out of memory in kp", CLEANUP);

    initialize(kp, env, sol, node);

    if (node->profit > sol->val)
        save_opt_sol(kp, env, sol, node);

    if (node->cap <= SOLVER_ZEROPLUS)
        goto CLEANUP;

    rval = dfs(kp, env, sol, node);
    check_rval(rval, "Branching failed in MT1", CLEANUP);

    if (env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("kp   : Number of branch nodes: %d\n", node->bbcount);
    }

CLEANUP:

    if (node->x)
        free(node->x);
    if (node)
        free(node);

    return rval;
}

static int
initialize(kp_prob *kp, kp_env *env, kp_sol *sol, branch_node *node)
{
    int rval = 0;
    int j;
    double min;

    kp_param *param = env->param;

    memset(node->x, -1, kp->n * sizeof(int));

    node->level      = -1;
    node->bbcount    = 0;
    node->branch     = 0;
    node->unbranch   = 0;
    node->profit     = 0.0;
    node->cap        = kp->cap * (1 + param->epsilon);
    node->upperbound = kp->upperbound;

    kp->sol_stat = SOLVER_FEAS;

    if (node->cap > 0.0)
    {
        param->p_epsilon = param->epsilon * sol->val;

        min                  = BIGDOUBLE;
        kp->min_w[kp->n - 1] = min;
        for (j = kp->n - 2; j >= 0; j--)
        {
            if (kp->w[j + 1] < min)
                min = kp->w[j + 1];
            kp->min_w[j] = min;
        }
    }

    return rval;
}

static int
dfs(kp_prob *kp, kp_env *env, kp_sol *sol, branch_node *node)
{
    int rval = 0;
    double oldbound;

    oldbound = node->upperbound;

    node->level++;
    node->bbcount++;

    // Branch to 1
    if (kp->w[node->level] <= node->cap)
    {
        node->branch = 1;
        rval         = branch(kp, env, sol, node);
        check_rval(rval, "branch failed\n", CLEANUP);

        if (node->cap == 0 || node->level == kp->n - 1)
        {
            if (node->profit > sol->val)
            {
                save_opt_sol(kp, env, sol, node);

                if (sol->val >= kp->upperbound)
                {
                    printf("Optimal sol found\n");
                    goto CLEANUP;
                }
            }
        }
        // else if (node->upperbound > sol->val)
        else if (node->upperbound > sol->val &&
                 node->profit + floor(node->cap * kp->p[node->level + 1] /
                                      kp->w[node->level + 1]) >
                 sol->val)
        {
            if (node->cap >= kp->min_w[node->level] &&
                node->bbcount <= env->max_counter)
            {
                rval = dfs(kp, env, sol, node);
                check_rval(rval, "dfs failed\n", CLEANUP);

                if (kp->sol_stat == SOLVER_OPT)
                    goto CLEANUP;
            }
        }
        node->unbranch = 1;
        rval           = unbranch(kp, env, node);

        node->upperbound = oldbound;
    }

    // Branch to 0
    if (node->level < kp->n - 1 && node->cap >= kp->min_w[node->level])
    {
        node->branch = -1;
        rval         = branch(kp, env, sol, node);
        check_rval(rval, "branch failed\n", CLEANUP);

        // if (node->upperbound > kp->sol->val)
        if (node->upperbound > sol->val &&
            node->profit +
            floor(node->cap * kp->p[node->level + 1] / kp->w[node->level + 1]) >
            sol->val)
        {
            if (node->cap >= kp->min_w[node->level] &&
                node->bbcount <= env->max_counter)
            {
                rval = dfs(kp, env, sol, node);
                check_rval(rval, "dfs failed\n", CLEANUP);

                if (kp->sol_stat == SOLVER_OPT)
                    goto CLEANUP;
            }
        }
        node->unbranch = -1;
        rval           = unbranch(kp, env, node);

        node->upperbound = oldbound;
    }

CLEANUP:

    node->level--;

    return rval;
}

int
branch(kp_prob *kp, kp_env *env, kp_sol *sol, branch_node *node)
{
    int rval = 0;

    if (node->branch > 0)
    {
        if (node->x[node->level] != -1)
            fprintf(stderr, "Error in branching\n");
        node->profit += kp->p[node->level];
        node->cap -= kp->w[node->level];
        node->x[node->level] = 1;

        if (sol->val < node->profit)
            save_opt_sol(kp, env, sol, node);
    }
    else if (node->branch < 0)
    {
        if (node->x[node->level] != -1)
            fprintf(stderr, "Error in branching\n");
        node->x[node->level] = 0;
    }
    else
    {
        fprintf(stderr, "Error in branching\n");
    }

    node->upperbound = get_upperbound(kp, env, node);

    node->branch = 0;

    return rval;
}

static int
unbranch(kp_prob *kp, kp_env *env, branch_node *node)
{
    int rval = 0;

    if (node->unbranch > 0)
    {
        if (node->x[node->level] != 1)
            fprintf(stderr, "Error in branching\n");
        node->cap += kp->w[node->level];
        node->profit -= kp->p[node->level];
        node->x[node->level] = -1;
    }
    else if (node->unbranch < 0)
    {
        if (node->x[node->level] != 0)
            fprintf(stderr, "Error in branching\n");
        node->x[node->level] = -1;
    }
    else
    {
        fprintf(stderr, "Error in branching\n");
    }

    node->unbranch = 0;

    return rval;
}

static int
save_opt_sol(kp_prob *kp, kp_env *env, kp_sol *sol, branch_node *node)
{
    int rval = 0;
    int i;

    sol->tot_n  = 0;
    sol->val    = 0.0;
    sol->weight = 0.0;
    sol->ns     = 0;
    memset(sol->sposition, kp->n, kp->n * sizeof(int));

    for (i = 0; i < kp->n; i++)
    {
        if (node->x[i] == 1)
        {
            sol->val += kp->p[i];
            sol->weight += kp->w[i];
            sol->selected[kp->ord_ind[i]] = 1;
            sol->ns++;
            ;
        }
        else if (!node->x[i])
        {
            sol->selected[kp->ord_ind[i]] = 0;
        }
        else if (node->x[i] == -1)
        {
            sol->selected[kp->ord_ind[i]] = 0;
        }
    }

    printf("Save %f\n", sol->val);

    if (sol->val >= kp->upperbound || kp->sol_stat == SOLVER_OPT)
        kp->sol_stat = SOLVER_OPT;
    else
        kp->sol_stat = SOLVER_FEAS;

    return rval;
}

static double
get_upperbound(kp_prob *kp, kp_env *env, branch_node *node)
{
    double upperbound = 0, tmpbound, tmpprofit;
    kp_param *param   = env->param;

    if (node->level < kp->n - 2)
    {
        if (env->p_integer)
            tmpprofit = ceil((kp->w[node->level + 1] - node->cap) *
                             kp->p[node->level] / kp->w[node->level]);
        else
            tmpprofit = (kp->w[node->level + 1] - node->cap) *
                        kp->p[node->level] / kp->w[node->level] +
                        param->epsilon;
        upperbound = node->profit + kp->p[node->level + 1] - tmpprofit;

        if (env->p_integer)
            tmpprofit =
            floor(node->cap * kp->p[node->level + 2] / kp->w[node->level + 2]);
        else
            tmpprofit =
            node->cap * kp->p[node->level + 2] / kp->w[node->level + 2] -
            param->epsilon;

        tmpbound = node->profit + tmpprofit;
        if (upperbound < tmpbound)
            upperbound = tmpbound;
    }
    else if (node->level == kp->n - 2)
    {
        if (env->p_integer)
            tmpprofit = ceil((kp->w[node->level + 1] - node->cap) *
                             kp->p[node->level] / kp->w[node->level]);
        else
            tmpprofit = (kp->w[node->level + 1] - node->cap) *
                        kp->p[node->level] / kp->w[node->level] +
                        param->epsilon;
        upperbound = node->profit + kp->p[node->level + 1] - tmpprofit;
    }
    else if (node->level == kp->n - 1)
    {
        upperbound = node->profit;
    }

    return upperbound;
}
