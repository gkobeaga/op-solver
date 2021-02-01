#include "op-solver.h"
#include "ip/exact/bac/bac.h"
#include "cp/exact/bac/bac.h"

int
ip_branch_dfs(ip_prob *ip, ip_exact_bac_env *bac_env, void *sol)
{
    int rval = 0;
    int prune;
    int nadded;
    int infeasible1 = 0;
    double bound_sense;

    double int_val;
    int sol_improved;

    ip_branch *branch             = NULL;
    ip_exact_bac_stats *bac_stats = bac_env->stats;
    ip_exact_bac_param *bac_param = bac_env->param;

    bac_env->depth++;

    if (ip->lowerboundG >= ip->upperboundG)
        goto CLEANUP;

    rval = stats_start(bac_stats->branch_node);
    check_rval(rval, "stats_start failed", CLEANUP);

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("\n");
        ip_print_branch_history(bac_env);
    }

    /**************************************************************************/
    /* Check if solution is valid                                             */
    /**************************************************************************/
    if (bac_env->check_sol(bac_env->lp->sol))
    {
        assert(bac_env->lp->sol->integral);
        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("Integral solution detected \n");

        int_val = bac_env->lp->sol->val;

        cp_cut *cuts = NULL;

        cp_sol *tmpsol = NULL;
        cp_get_sol_from_graph(bac_env->orig_prob, bac_env->lp->sol->graph,
                              &tmpsol);
        cuts = cp_conv_cut_sol2connect(tmpsol);
        rval =
        cp_add_lp_cuts(bac_env->orig_prob, bac_env->orig_env, &cuts, &nadded);

        if (bac_env->lp->status == SOLVER_LP_INFEASIBLE)
        {
            rval = stats_start(bac_stats->verify_branch_infeas);
            check_rval(rval, "stats_start failed", CLEANUP);
            rval = bac_env->verify_infeas(bac_env->orig_prob, bac_env->orig_env,
                                          &prune);
            check_rval(rval, "ip_verify_infeasible failed", CLEANUP);
            rval = stats_stop(bac_stats->verify_branch_infeas, prune);
            check_rval(rval, "stats_stop failed", CLEANUP);
            if (prune)
            {
                rval = stats_stop(bac_stats->branch_node, 1);
                check_rval(rval, "stats_stop failed", CLEANUP);
                if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
                    printf("PRUNE SEARCH - infeasible LP\n");
                ip_write_exact_bac_stats(bac_env);
                goto CLEANUP;
            }
            else
            {
                rval = stats_stop(bac_stats->branch_node, 0);
                check_rval(rval, "stats_stop failed", CLEANUP);
                check_assert(prune == 1,
                             "exact pricing did not verify an infeasible LP",
                             CLEANUP);
            }
        }

        sol_improved = 0;
        if (ip->sense == SOLVER_OPT_SENSE_MIN)
        {
            if (int_val < ip->upperboundG - SOLVER_ZEROPLUS)
            {
                bac_env->upperboundN = int_val;
                ip->upperboundG      = int_val;

                sol_improved = 1;
                cp_copy_sol(tmpsol, sol);
            }
        }
        else
        {
            if (int_val > ip->lowerboundG + SOLVER_ZEROPLUS)
            {
                bac_env->lowerboundN = int_val;
                ip->lowerboundG      = int_val;
                sol_improved         = 1;
                cp_copy_sol(tmpsol, sol);
            }
        }

        if (ip->lowerboundG >= ip->upperboundG)
            goto CLEANUP;
    }

    /**************************************************************************/
    /* Pricing Loop                                                           */
    /**************************************************************************/
    rval = stats_start(bac_stats->pricing_loop);
    check_rval(rval, "stats_start failed", CLEANUP);
    rval =
    bac_env->pricing_loop(bac_env->orig_prob, bac_env->orig_env, &nadded);
    check_rval(rval, "ip_pricing_loop failed", CLEANUP);
    rval = stats_stop(bac_stats->pricing_loop, nadded);
    check_rval(rval, "stats_stop failed", CLEANUP);

    if (bac_env->lowerboundN >= bac_env->upperboundN - bac_param->pruning_tol)
    {
        rval = stats_start(bac_stats->verify_branch_prune);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval =
        bac_env->verify_prune(bac_env->orig_prob, bac_env->orig_env, &prune);
        check_rval(rval, "ip_verify_prune failed", CLEANUP);
        ip_write_exact_bac_stats(bac_env);
        if (prune)
        {
            rval = stats_stop(bac_stats->verify_branch_prune, 1);
            check_rval(rval, "stats_stop failed", CLEANUP);
            rval = stats_stop(bac_stats->branch_node, 1);
            check_rval(rval, "stats_stop failed", CLEANUP);
            rval = 0;
            goto CLEANUP;
        }
        else
        {
            rval = stats_stop(bac_stats->verify_branch_prune, 0);
            check_rval(rval, "stats_stop failed", CLEANUP);
            rval = stats_stop(bac_stats->branch_node, 0);
            check_rval(rval, "stats_stop failed", CLEANUP);
            check_assert(prune == 1,
                         "exact pricing did not verify an infeasible LP",
                         CLEANUP);
        }
    }

    /**************************************************************************/
    /* Separation Loop                                                        */
    /**************************************************************************/
    rval = stats_start(bac_stats->sep_loop);
    check_rval(rval, "stats_start failed", CLEANUP);
    rval = bac_env->sep_loop(bac_env->orig_prob, bac_env->orig_env, &nadded);
    check_rval(rval, "failed", CLEANUP);
    rval = stats_stop(bac_stats->sep_loop, nadded);
    check_rval(rval, "stats_stop failed", CLEANUP);

    if (ip->lowerboundG >= ip->upperboundG)
        goto CLEANUP;

    if (bac_env->lp->status == SOLVER_LP_INFEASIBLE)
    {
        rval = stats_start(bac_stats->verify_branch_infeas);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval =
        bac_env->verify_infeas(bac_env->orig_prob, bac_env->orig_env, &prune);
        check_rval(rval, "ip_verify_infeasible failed", CLEANUP);
        rval = stats_stop(bac_stats->verify_branch_infeas, prune);
        check_rval(rval, "stats_stop failed", CLEANUP);
        if (prune)
        {
            rval = stats_stop(bac_stats->branch_node, 1);
            check_rval(rval, "stats_stop failed", CLEANUP);
            if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
                printf("PRUNE SEARCH - infeasible LP\n");
            ip_write_exact_bac_stats(bac_env);
            goto CLEANUP;
        }
        else
        {
            rval = stats_stop(bac_stats->branch_node, 0);
            check_rval(rval, "stats_stop failed", CLEANUP);
            check_assert(prune == 1,
                         "exact pricing did not verify an infeasible LP",
                         CLEANUP);
        }
    }
    else if (stats_get_current_time(bac_stats->total) > bac_param->time_limit)
    {
        rval = stats_stop(bac_stats->branch_node, 1);
        check_rval(rval, "stats_stop failed", CLEANUP);
        rval = 0;
        goto CLEANUP;
    }

    if (ip->sense == SOLVER_OPT_SENSE_MIN)
    {
        if (ip->upperboundG > bac_env->upperboundN - SOLVER_ZEROPLUS)
            ip->upperboundG = bac_env->upperboundN;
        if (bac_env->upperboundN > ip->upperboundG - SOLVER_ZEROPLUS)
            bac_env->upperboundN = ip->upperboundG;
    }
    else
    {
        if (ip->lowerboundG < bac_env->lowerboundN - SOLVER_ZEROPLUS)
            ip->lowerboundG = bac_env->lowerboundN;
        if (bac_env->lowerboundN < ip->lowerboundG - SOLVER_ZEROPLUS)
            bac_env->lowerboundN = ip->lowerboundG;
    }

    if (bac_env->lowerboundN >= bac_env->upperboundN - bac_param->pruning_tol)
    {
        rval = stats_start(bac_stats->verify_branch_prune);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval =
        bac_env->verify_prune(bac_env->orig_prob, bac_env->orig_env, &prune);
        check_rval(rval, "ip_verify_prune failed", CLEANUP);
        ip_write_exact_bac_stats(bac_env);
        if (prune)
        {
            rval = stats_stop(bac_stats->verify_branch_prune, 1);
            check_rval(rval, "stats_stop failed", CLEANUP);
            rval = stats_stop(bac_stats->branch_node, 1);
            check_rval(rval, "stats_stop failed", CLEANUP);
            rval = 0;
            goto CLEANUP;
        }
        else
        {
            rval = stats_stop(bac_stats->verify_branch_prune, 0);
            check_rval(rval, "stats_stop failed", CLEANUP);
            rval = stats_stop(bac_stats->branch_node, 0);
            check_rval(rval, "stats_stop failed", CLEANUP);
            check_assert(prune == 1,
                         "exact pricing did not verify an infeasible LP",
                         CLEANUP);
        }
    }

    /**************************************************************************/
    /* Obtain Heuristic Solution form LP Value                                */
    /**************************************************************************/
    rval = bac_env->xheur(bac_env->orig_prob, bac_env->orig_env, sol, &nadded);
    check_rval(rval, " failed", CLEANUP);

    if (ip->lowerboundG >= ip->upperboundG)
        goto CLEANUP;

    if (bac_env->lp->status == SOLVER_LP_INFEASIBLE)
    {
        rval = stats_start(bac_stats->verify_branch_infeas);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval =
        bac_env->verify_infeas(bac_env->orig_prob, bac_env->orig_env, &prune);
        check_rval(rval, "ip_verify_infeasible failed", CLEANUP);
        rval = stats_stop(bac_stats->verify_branch_infeas, prune);
        check_rval(rval, "stats_stop failed", CLEANUP);
        if (prune)
        {
            rval = stats_stop(bac_stats->branch_node, 1);
            check_rval(rval, "stats_stop failed", CLEANUP);
            if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
                printf("PRUNE SEARCH - infeasible LP\n");
            ip_write_exact_bac_stats(bac_env);
            goto CLEANUP;
        }
        else
        {
            rval = stats_stop(bac_stats->branch_node, 0);
            check_rval(rval, "stats_stop failed", CLEANUP);
            check_assert(prune == 1,
                         "exact pricing did not verify an infeasible LP",
                         CLEANUP);
        }
    }

    if (bac_env->lowerboundN >= bac_env->upperboundN - bac_param->pruning_tol)
    {
        rval = stats_start(bac_stats->verify_branch_prune);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval =
        bac_env->verify_prune(bac_env->orig_prob, bac_env->orig_env, &prune);
        check_rval(rval, "ip_verify_prune failed", CLEANUP);
        ip_write_exact_bac_stats(bac_env);
        if (prune)
        {
            rval = stats_stop(bac_stats->verify_branch_prune, 1);
            check_rval(rval, "stats_stop failed", CLEANUP);
            rval = stats_stop(bac_stats->branch_node, 1);
            check_rval(rval, "stats_stop failed", CLEANUP);
            rval = 0;
            goto CLEANUP;
        }
        else
        {
            rval = stats_stop(bac_stats->verify_branch_prune, 0);
            check_rval(rval, "stats_stop failed", CLEANUP);
            rval = stats_stop(bac_stats->branch_node, 0);
            check_rval(rval, "stats_stop failed", CLEANUP);
            check_assert(prune == 1,
                         "exact pricing did not verify an infeasible LP",
                         CLEANUP);
        }
    }

    /**************************************************************************/
    /* Check if Solution is Integral                                          */
    /**************************************************************************/

    if (bac_env->check_sol(bac_env->lp->sol))
    {
        assert(bac_env->lp->sol->integral);
        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("Integral solution detected \n");
        // if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)

        int_val = bac_env->lp->sol->val;

        cp_cut *cuts = NULL;

        cp_sol *tmpsol = NULL;
        cp_get_sol_from_graph(bac_env->orig_prob, bac_env->lp->sol->graph,
                              &tmpsol);
        cuts = cp_conv_cut_sol2connect(tmpsol);
        rval =
        cp_add_lp_cuts(bac_env->orig_prob, bac_env->orig_env, &cuts, &nadded);

        if (nadded && ip->lp->status == SOLVER_LP_SUCCESS &&
            (bac_env->upperboundN <
             bac_env->lowerboundN + bac_env->param->pruning_tol))
        {
            int edge_added;
            rval = cp_add_badvars(bac_env->orig_prob, bac_env, &edge_added);
        }

        sol_improved = 0;
        if (ip->sense == SOLVER_OPT_SENSE_MIN)
        {
            if (int_val < ip->upperboundG - SOLVER_ZEROPLUS)
            {
                bac_env->upperboundN = int_val;
                ip->upperboundG      = int_val;

                sol_improved = 1;
                cp_copy_sol(tmpsol, sol);
            }
        }
        else
        {
            if (int_val > ip->lowerboundG + SOLVER_ZEROPLUS)
            {
                bac_env->lowerboundN = int_val;
                ip->lowerboundG      = int_val;
                sol_improved         = 1;
                cp_copy_sol(tmpsol, sol);
            }
        }

        if (ip->lowerboundG >= ip->upperboundG)
            goto CLEANUP;

        if ((sol_improved || nadded) &&
            bac_env->lowerboundN >=
            bac_env->upperboundN - bac_param->pruning_tol)
        {
            rval = stats_start(bac_stats->pricing_loop);
            check_rval(rval, "stats_start failed", CLEANUP);
            rval = bac_env->pricing_loop(bac_env->orig_prob, bac_env->orig_env,
                                         &nadded);
            check_rval(rval, "ip_pricing_loop failed", CLEANUP);
            rval = stats_stop(bac_stats->pricing_loop, nadded);
            check_rval(rval, "stats_stop failed", CLEANUP);
        }

        if (bac_env->lowerboundN >=
            bac_env->upperboundN - bac_param->pruning_tol)
        {
            rval = bac_env->verify_prune(bac_env->orig_prob, bac_env->orig_env,
                                         &prune);
            check_rval(rval, "ip_verify_prune failed", CLEANUP);
            ip_write_exact_bac_stats(bac_env);
            if (prune)
            {
                rval = stats_stop(bac_stats->branch_node, 1);
                check_rval(rval, "stats_start failed", CLEANUP);
                printf("with new tour, the node can be pruned\n");
                ip_write_exact_bac_stats(bac_env);
                rval = 0;
                goto CLEANUP;
            }
            else
            {
                rval = stats_stop(bac_stats->branch_node, 0);
                check_rval(rval, "stats_start failed", CLEANUP);
                check_assert(prune == 1, "could not verify the pruning",
                             CLEANUP);
            }
        }
    }

    rval = stats_stop(bac_stats->branch_node, 1);
    check_rval(rval, "stats_start failed", CLEANUP);

    /**************************************************************************/
    /* Update bounds and write stats                                          */
    /**************************************************************************/

    if (ip->sense == SOLVER_OPT_SENSE_MIN)
    {
        if (ip->upperboundG > bac_env->upperboundN + SOLVER_ZEROPLUS)
            ip->upperboundG = bac_env->upperboundN;
        if (bac_env->upperboundN > ip->upperboundG + SOLVER_ZEROPLUS)
            bac_env->upperboundN = ip->upperboundG;
        bound_sense = bac_env->lowerboundN;
    }
    else
    {
        if (ip->lowerboundG < bac_env->lowerboundN - SOLVER_ZEROPLUS)
            ip->lowerboundG = bac_env->lowerboundN;
        if (bac_env->lowerboundN < ip->lowerboundG - SOLVER_ZEROPLUS)
            bac_env->lowerboundN = ip->lowerboundG;
        bound_sense = bac_env->upperboundN;
    }

    if (bac_env->depth == 1)
    {
        if (ip->sense == SOLVER_OPT_SENSE_MIN &&
            ip->lowerboundG < bac_env->lowerboundN - SOLVER_ZEROPLUS)
            ip->lowerboundG = bac_env->lowerboundN;
        else if (ip->sense == SOLVER_OPT_SENSE_MAX &&
                 ip->upperboundG > bac_env->upperboundN + SOLVER_ZEROPLUS)
            ip->upperboundG = bac_env->upperboundN;
    }

    ip_write_exact_bac_stats(bac_env);

    /**************************************************************************/
    /* Find Edge to Branch                                                    */
    /**************************************************************************/

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        printf("Find branch object ...\n");

    rval = stats_start(bac_stats->find_branch);
    check_rval(rval, "stats_start failed", CLEANUP);
    rval = bac_env->find_branch(bac_env->orig_prob, bac_env->orig_env, &branch);
    check_rval(rval, "failed", CLEANUP);
    rval = stats_stop(bac_stats->find_branch, branch->edge != NULL);
    check_rval(rval, "stats_start failed", CLEANUP);

    /**************************************************************************/
    /* Up-Side Branch                                                         */
    /**************************************************************************/
    branch->rhs = 1;
    // branch->rhs = rng_bernoulli(0.5);

    rval = stats_start(bac_stats->exec_branch);
    check_rval(rval, "stats_start failed", CLEANUP);
    rval = ip_exec_branch(ip, bac_env, branch);
    check_rval(rval, "failed", CLEANUP);
    rval = stats_stop(bac_stats->exec_branch, 1);
    check_rval(rval, "stats_stop failed", CLEANUP);
    if (bac_env->lp->status == SOLVER_LP_INFEASIBLE)
    {
        rval = stats_start(bac_stats->verify_branch_infeas);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval =
        bac_env->verify_infeas(bac_env->orig_prob, bac_env->orig_env, &prune);
        check_rval(rval, "ip_verify_infeasible failed", CLEANUP);
        rval = stats_stop(bac_stats->verify_branch_infeas, prune);
        check_rval(rval, "stats_stop failed", CLEANUP);
        if (prune)
        {
            if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
                printf("PRUNE SIDE - infeasible LP\n");
            if (bac_env->depth == 1)
                infeasible1 = 1;
            rval = 0;
        }
        else
        {
            check_assert(prune == 1, "could not verify the pruning", CLEANUP);
        }
    }
    else if (bac_env->lp->status == SOLVER_LP_SUCCESS)
    {
        rval = ip_branch_dfs(ip, bac_env, sol);
        check_rval(rval, " failed", CLEANUP);
    }
    else
    {
        fprintf(stderr, "error in DFS\n");
        rval = 1;
        goto CLEANUP;
    }

    rval = stats_start(bac_stats->exec_unbranch);
    check_rval(rval, "stats_start failed", CLEANUP);
    rval = ip_exec_unbranch(ip, bac_env);
    check_rval(rval, "failed", CLEANUP);
    rval = stats_stop(bac_stats->exec_unbranch, 1);
    check_rval(rval, "stats_stop failed", CLEANUP);

    if (ip->lowerboundG >= ip->upperboundG)
        goto CLEANUP;

    if (bac_env->lp->status == SOLVER_LP_INFEASIBLE)
    {
        rval = stats_start(bac_stats->verify_branch_infeas);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval =
        bac_env->verify_infeas(bac_env->orig_prob, bac_env->orig_env, &prune);
        check_rval(rval, "ip_verify_infeasible failed", CLEANUP);
        rval = stats_stop(bac_stats->verify_branch_infeas, prune);
        check_rval(rval, "stats_stop failed", CLEANUP);
        if (prune)
        {
            // if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("PRUNE BOTH SIDES - infeasible unbranched LP\n");
            rval = 0;
            if (ip->sense == SOLVER_OPT_SENSE_MIN)
                bac_env->lowerboundN = bound_sense;
            else
                bac_env->upperboundN = bound_sense;
            ip_free_branch(&branch);
            goto CLEANUP;
        }
        else
        {
            fprintf(stderr, "exact pricing did not verify an infeasible LP\n");
            rval = 1;
            goto CLEANUP;
        }
    }

    if (ip->sense == SOLVER_OPT_SENSE_MIN)
        bac_env->lowerboundN = bound_sense;
    else
        bac_env->upperboundN = bound_sense;

    if (stats_get_total_time(bac_stats->total) > bac_param->time_limit)
    {
        rval = 0;
        goto CLEANUP;
    }

    /**************************************************************************/
    /* Down-Side Branch                                                       */
    /**************************************************************************/
    branch->rhs = !branch->rhs;

    rval = stats_start(bac_stats->exec_branch);
    check_rval(rval, "stats_start failed", CLEANUP);
    rval = ip_exec_branch(ip, bac_env, branch);
    rval = stats_stop(bac_stats->exec_branch, 1);
    check_rval(rval, "stats_stop failed", CLEANUP);
    if (bac_env->lp->status == SOLVER_LP_INFEASIBLE)
    {
        fprintf(stderr, "branched lp was infeasible\n");
        rval = stats_start(bac_stats->verify_branch_infeas);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval =
        bac_env->verify_infeas(bac_env->orig_prob, bac_env->orig_env, &prune);
        check_rval(rval, "ip_verify_infeasible failed", CLEANUP);
        rval = stats_stop(bac_stats->verify_branch_infeas, prune);
        check_rval(rval, "stats_stop failed", CLEANUP);
        if (prune)
        {
            if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
                printf("PRUNE SIDE - infeasible LP\n");
            rval = 0;
        }
        else
        {
            check_assert(prune == 1, "could not verify the pruning", CLEANUP);
        }
    }
    else if (bac_env->lp->status == SOLVER_LP_SUCCESS)
    {
        rval = ip_branch_dfs(ip, bac_env, sol);
        check_rval(rval, "failed", CLEANUP);
    }
    else
    {
        fprintf(stderr, "error in DFS\n");
        rval = 1;
        goto CLEANUP;
    }

    rval = stats_start(bac_stats->exec_unbranch);
    check_rval(rval, "stats_start failed", CLEANUP);
    rval = ip_exec_unbranch(ip, bac_env);
    check_rval(rval, "failed", CLEANUP);
    rval = stats_stop(bac_stats->exec_unbranch, 1);
    check_rval(rval, "stats_stop failed", CLEANUP);

    if (ip->lowerboundG >= ip->upperboundG)
        goto CLEANUP;

    if (bac_env->lp->status == SOLVER_LP_INFEASIBLE)
    {
        fprintf(stderr, "unbranched lp was infeasible\n");
        rval = stats_start(bac_stats->verify_branch_infeas);
        check_rval(rval, "stats_start failed", CLEANUP);
        rval =
        bac_env->verify_infeas(bac_env->orig_prob, bac_env->orig_env, &prune);
        check_rval(rval, "ip_verify_infeasible failed", CLEANUP);
        rval = stats_stop(bac_stats->verify_branch_infeas, prune);
        check_rval(rval, "stats_stop failed", CLEANUP);
        if (prune)
        {
            printf("NOTE - infeasible unbranched LP\n");
            if (bac_env->depth == 1 && infeasible1)
                ip->infeasible = 1;
            goto CLEANUP;
        }
        else
        {
            check_assert(prune == 1, "could not verify the pruning", CLEANUP);
        }
    }

    if (ip->sense == SOLVER_OPT_SENSE_MIN)
        bac_env->lowerboundN = bound_sense;
    else
        bac_env->upperboundN = bound_sense;

    if (rval == 0 && bac_env->depth == 1 &&
        stats_get_total_time(bac_stats->total) <= bac_param->time_limit)
    {
        if (ip->sense == SOLVER_OPT_SENSE_MIN)
            ip->lowerboundG = ip->upperboundG;
        else
            ip->upperboundG = ip->lowerboundG;
    }

CLEANUP:
    if (branch)
        ip_free_branch(&branch);
    bac_env->depth--;

    return rval;
}
