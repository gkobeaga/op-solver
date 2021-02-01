#include "op-solver.h"
#include "ip/ip.h"
#include "ip/exact/bac/bac.h"

typedef struct sbranch_item
{
    int i;
    int lpind;
    double val;
} sbranch_item;

ip_branch *
ip_create_branch(void)
{
    ip_branch *branch = malloc(sizeof(ip_branch));
    branch->depth     = 0;
    branch->rhs       = 0;
    branch->edge      = NULL;
    branch->sense     = 'X';
    branch->lprow     = 0;
    return branch;
}

void
ip_free_branch(ip_branch **branch)
{
    if (*branch)
    {
        (*branch)->depth = 0;
        (*branch)->rhs   = 0;
        (*branch)->edge  = NULL;
        (*branch)->sense = 'X';
        free(*branch);
        *branch = NULL;
    }
}

int
ip_exec_branch(ip_prob *ip, ip_exact_bac_env *bac_env, ip_branch *branch)
{
    int rval = 0;
    graph_arc *edge;
    lp_data *data = NULL;

    ip_exact_bac_stats *bac_stats = bac_env->stats;
    lp_prob *lp                   = bac_env->lp;

    check_null(branch, "ip_exec_branch called without a ip_branch failed\n",
               CLEANUP);

    edge = branch->edge;
    check_null(edge, "ip_exec_branch has invalid edge", CLEANUP);
    check_rval(edge->fixed, "branching edge is fixed to 1 in the LP", CLEANUP);

    rval = (int)(edge->branch);
    check_rval(rval, "branching edge has already been branched", CLEANUP);

    if (branch->rhs)
    {
        rval = lp_set_bnd(lp, lp->graph->nv + edge->ind, 'L', 1.0);
        check_rval(rval, "lp_set_bnd failed", CLEANUP);
        edge->branch = bac_env->depth + 1;
        assert(edge->branch);
    }
    else
    {
        rval = lp_set_bnd(lp, lp->graph->nv + edge->ind, 'U', 0.0);
        check_rval(rval, "lp_set_bnd failed", CLEANUP);
        edge->branch = -(bac_env->depth + 1);
    }

    rval = lp_opt_dual(lp);

    if (ip->lowerboundG >= ip->upperboundG)
        goto CLEANUP;

    if (lp->status == SOLVER_LP_INFEASIBLE)
    {
        rval = stats_start(bac_stats->recover_infeas);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval = bac_env->recover_infeas(bac_env->orig_prob, bac_env->orig_env);
        check_rval(rval, " failed", CLEANUP);
        if (lp->status == SOLVER_LP_INFEASIBLE)
        {
            rval = stats_stop(bac_stats->recover_infeas, 0);
            check_rval(rval, "failed", CLEANUP);
            if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
            {
                if (lp->status == SOLVER_LP_INFEASIBLE)
                    printf("BRANCH LP INFEASIBLE\n");
            }
            goto CLEANUP;
        }
        else

        {
            rval = stats_stop(bac_stats->recover_infeas, 1);
            check_rval(rval, "failed", CLEANUP);
        }
    }
    else if (rval)
    {
        fprintf(stderr, " failed\n");
        rval = 1;
        goto CLEANUP;
    }

    rval = lp_update_sol(lp, lp->sol);
    check_rval(rval, " failed", CLEANUP);

CLEANUP:

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf(
        "(%d/%d) :   Branch edge (%d, %d) to   %d : LB %.0f : LP %.2f: UB "
        "%.0f : GAP %.2f : time %ld\n",
        bac_env->stats->exec_branch->count_active -
        bac_env->stats->exec_unbranch->count_active + 1,
        bac_env->stats->exec_branch->count_active + 1, branch->edge->tail->i,
        branch->edge->head->i, branch->rhs, ip->lowerboundG,
        lp->status == SOLVER_LP_INFEASIBLE ? 0 : bac_env->lp->sol->val,
        ip->upperboundG,

        (ip->upperboundG - ip->lowerboundG) / ip->upperboundG * 100,
        stats_get_total_time(bac_env->stats->total));
    }

    if (rval == 0 || lp->status == SOLVER_LP_INFEASIBLE)
    {
        if (bac_env->history_space < bac_env->history_depth + 1)
            realloc_scale(bac_env->history, bac_env->history_space,
                          bac_env->history_depth + 1, 1.3);
        check_null(bac_env->history, "realloc failed", CLEANUP);
        bac_env->history[bac_env->history_depth] = ip_create_branch();
        bac_env->history[bac_env->history_depth]->depth =
        bac_env->history_depth + 1;
        bac_env->history[bac_env->history_depth]->edge  = branch->edge;
        bac_env->history[bac_env->history_depth]->rhs   = branch->rhs;
        bac_env->history[bac_env->history_depth]->sense = branch->sense;
        bac_env->history[bac_env->history_depth]->lprow = branch->lprow;
        bac_env->history_depth++;
    }
    lp_free_data(&data);
    return rval;
}

int
ip_exec_unbranch(ip_prob *ip, ip_exact_bac_env *bac_env)
{
    int rval = 0;

    ip_exact_bac_stats *bac_stats = bac_env->stats;
    lp_prob *lp                   = bac_env->lp;

    int depth         = bac_env->depth;
    ip_branch *branch = NULL;
    graph_arc *edge;

    rval = (int)(depth <= 0);
    check_rval(rval, "called at depth 0", CLEANUP);
    rval = (int)(bac_env->history[depth - 1]->depth != depth);
    check_rval(rval, "branchhistory is corrupted", CLEANUP);

    branch = bac_env->history[depth - 1];

    edge = branch->edge;

    check_null(edge, "ERROR: unbranching 1-edge is not in iP", CLEANUP);

    if (branch->rhs)
    {
        if (edge->branch <= 0)
        {
            rval = 1;
            check_rval(rval, "unbranching 1-edge not branched to 1", CLEANUP);
        }

        rval = lp_set_bnd(lp, lp->graph->nv + edge->ind, 'L', 0.0);
        check_rval(rval, "failed", CLEANUP);
    }
    else
    {
        if (edge->branch >= 0)
        {
            rval = 1;
            check_rval(rval, "unbranching 0-edge not branched to 0", CLEANUP);
        }

        rval = lp_set_bnd(lp, lp->graph->nv + edge->ind, 'U', 1.0);
        check_rval(rval, "failed", CLEANUP);
    }
    edge->branch = 0;

    if (ip->lowerboundG >= ip->upperboundG)
        goto CLEANUP;

    rval = lp_opt_dual(lp);

    if (lp->status == SOLVER_LP_INFEASIBLE)
    {
        rval = stats_start(bac_stats->recover_infeas);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval = bac_env->recover_infeas(bac_env->orig_prob, bac_env->orig_env);
        check_rval(rval, "ip_recover_infeas failed", CLEANUP);
        if (lp->status == SOLVER_LP_INFEASIBLE)
        {
            rval = stats_stop(bac_stats->recover_infeas, 0);
            check_rval(rval, "stats_stop failed", CLEANUP);
            if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
            {
                printf("Problem is really infeasible\n");
            }
            goto CLEANUP;
        }
        else
        {
            rval = stats_stop(bac_stats->recover_infeas, 1);
            check_rval(rval, "stats_stop failed", CLEANUP);
        }
    }
    else if (rval)
    {
        fprintf(stderr, "failed\n");
        rval = 1;
        goto CLEANUP;
    }

    rval = lp_update_sol(lp, lp->sol);
    check_rval(rval, "failed", CLEANUP);

CLEANUP:

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf(
        "(%d/%d) : Unbranch edge (%d, %d) from %d : LB %.0f : LP %.2f: UB "
        "%.0f : GAP %.2f : time %ld\n",
        bac_env->stats->exec_branch->count_active -
        bac_env->stats->exec_unbranch->count_active + 1,
        bac_env->stats->exec_branch->count_active + 1, branch->edge->tail->i,
        branch->edge->head->i, branch->rhs, ip->lowerboundG,
        lp->status == SOLVER_LP_INFEASIBLE ? 0 : bac_env->lp->sol->val,
        ip->upperboundG,
        (ip->upperboundG - ip->lowerboundG) / ip->upperboundG * 100,
        stats_get_total_time(bac_env->stats->total));
    }

    if (rval == 0 || lp->status == SOLVER_LP_INFEASIBLE)
    {
        ip_free_branch(&(bac_env->history[bac_env->history_depth - 1]));
        bac_env->history_depth--;
    }
    return rval;
}
