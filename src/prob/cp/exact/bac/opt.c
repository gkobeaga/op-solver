#include "op-solver.h"

int
cp_opt_exact_bac(cp_prob *cp, cp_exact_bac_env *bac_env, cp_sol *sol)
{
    int rval = 0;

    cp_exact_bac_stats *bac_stats = bac_env->stats;

    int nadded;
    cp_sol *tmpsol = NULL;

    if (!sol->val)
    {
        bac_env->heur->param->init = SOLVER_CP_INIT_BEST3;
        rval                       = cp_init_sol(cp, bac_env->heur, sol);
        bac_env->heur->param->init = SOLVER_CP_INIT_RAND;

        /// First, we obtain an approximation of the TSP value of the data
        data_map *cpmap = cp->data->map;
        int *selected   = malloc(cp->n * sizeof(int));
        memset(selected, 1, cp->n * sizeof(int));
        data_emb_map(cp->data, selected);
        free(selected);

        tsp_prob *tsp         = tsp_create_prob(cp->data);
        tsp_init_env *tspinit = tsp_create_init_env();
        tsp_init_sol(tsp, tspinit, tsp->sol);
        tsp_free_init_env(&tspinit);

        /// Get the edge set in NN-10. Used in the Lin-Kernighan
        rval = data_get_k_nearest(tsp->data, 10);
        check_rval(rval, "", CLEANUP);

        rval = tsp_opt_heur(tsp, bac_env->heur->tsp, tsp->sol);
        check_rval(rval, "", CLEANUP);

        /// Define the probability to include a vertex in the initial solution.
        /// p = sqrt(d_0/TSP)
        bac_env->heur->param->pinit = sqrt(cp->data->cap / tsp->sol->val);
        data_free_map(&cp->data->map);
        cp->data->map = cpmap;
        tsp_free_prob(&tsp);

        /// Create the population of feasible solutions
        cp_pop *pop = cp_create_pop(cp, bac_env->heur->ea->param->pop_size);
        check_null(pop, "", CLEANUP);
        rval = cp_init_pop(cp, bac_env->heur, pop);
        check_rval(rval, "", CLEANUP);

        /// Use the EA4OP metaheuristic to find an initial solution
        rval = cp_opt_heur_ea(cp, bac_env->heur, pop, sol);
        check_rval(rval, "", CLEANUP);

        cp_free_pop(&pop);
    }

    rval = cp_init_exact_bac(cp, bac_env);
    check_rval(rval, "cp_init_exact_bac failed", CLEANUP);

    rval = stats_start(bac_stats->total);
    check_rval(rval, "stats_start failed", CLEANUP);
    rval = stats_start(bac_env->ip->stats->total);
    check_rval(rval, "stats_start failed", CLEANUP);

    lp_prob *lp = bac_env->ip->lp;

    if (bac_env->ip->lp->status == SOLVER_LP_INFEASIBLE)
    {
        printf("ip_init_prob reports infeasible LP subproblem\n");

        rval = bac_env->ip->recover_infeas(cp, bac_env);
        check_rval(rval, "recover_infeas failed", CLEANUP);
        if (bac_env->ip->lp->status != SOLVER_LP_INFEASIBLE)
        {
            printf("Couldn't verify infeasible IP\n");
            rval = 1;
            goto CLEANUP;
        }
        printf("ip_init_prob reports infeasible CP problem\n");
        goto DONE;
    }
    else if (rval)
    {
        fprintf(stderr, "Infeasible CP problem\n");
        goto CLEANUP;
    }

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("Inital LP has %d rows, %d columns, %d nonzeros. Val %f\n",
               lp_get_nrows(lp), lp_get_ncols(lp), lp_get_nnonzeros(lp),
               lp_get_objval(lp));
    }

    rval = bac_env->ip->sep_loop(cp, bac_env, &nadded);
    check_rval(rval, " failed", CLEANUP);

    if (bac_env->ip->lp->status == SOLVER_LP_INFEASIBLE)
    {
        printf("sep_locp reports an infeasible LP subproblem\n");

        rval = bac_env->ip->recover_infeas(cp, bac_env);
        check_rval(rval, "recover_infeas failed", CLEANUP);
        if (bac_env->ip->lp->status != SOLVER_LP_INFEASIBLE)
        {
            printf("Couldn't verify infeasible LP\n");
            rval = 1;
            goto CLEANUP;
        }

        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("Final LP has %d rows, %d columns, %d nonzeros\n",
                   lp_get_nrows(lp), lp_get_ncols(lp), lp_get_nnonzeros(lp));
        }

        goto DONE;
    }

    if (bac_env->ip->lowerboundN <
        bac_env->ip->upperboundN - bac_env->ip->param->pruning_tol)
    {
        tmpsol = cp_create_sol(cp);
        check_null(tmpsol, "cp_create_sol failed", CLEANUP);

        rval = stats_start(bac_stats->xheur_branch);
        check_rval(rval, "stats_start failed", CLEANUP);

        rval = bac_env->ip->xheur(cp, bac_env, tmpsol, &nadded);
        check_rval(rval, "ip->xheur failed", CLEANUP);

        rval = stats_stop(bac_stats->xheur_branch, nadded);
        check_rval(rval, "failed", CLEANUP);

        if (tmpsol->val > sol->val + SOLVER_ZEROPLUS)
        {
            if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
                printf("New lowerbound from primal heuristic: %.2f\n",
                       tmpsol->val);
            cp->sol_status = SOLVER_FEAS;
            cp_copy_sol(tmpsol, sol);
        }

        cp_free_sol(&tmpsol);
    }

    mpf_t lowbound, upbound, one;
    mpf_init(one);
    mpf_set_d(one, 1 - SOLVER_PRICE_MAXPENALTY);
    mpf_init(lowbound);
    mpf_set_d(lowbound, cp->ip->lowerboundG);
    mpf_init(upbound);

    rval = bac_env->ip->dual_bound(cp, bac_env, upbound);
    check_rval(rval, "ip->dual_bound failed", CLEANUP);

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("Dual upper bound: %f\n", mpf_get_d(upbound));
        printf("Primal upper bound: %f\n", bac_env->ip->lp->sol->val);
        printf("Primal lowerbound bound: %f\n", mpf_get_d(lowbound));
    }

    mpf_sub(upbound, upbound, one);
    mpf_clear(one);
    if (mpf_cmp(lowbound, upbound) > 0)
    {

        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("Final LP has %d rows, %d columns, %d nonzeros\n",
                   lp_get_nrows(lp), lp_get_ncols(lp), lp_get_nnonzeros(lp));
        }

        ip_write_exact_bac_stats(bac_env->ip);
        mpf_clear(lowbound);
        mpf_clear(upbound);
        goto DONE;
    }
    mpf_clear(lowbound);
    mpf_clear(upbound);

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("Final LP has %d rows, %d columns, %d nonzeros\n",
               lp_get_nrows(lp), lp_get_ncols(lp), lp_get_nnonzeros(lp));
    }

    ip_write_exact_bac_stats(bac_env->ip);
    rval = ip_branch_dfs(cp->ip, bac_env->ip, sol);
    check_rval(rval, "ip_branch_dfs failed", CLEANUP);
    ip_write_exact_bac_stats(bac_env->ip);

DONE:

    cp_write_exact_bac_stats(cp, bac_env);

    if (bac_env->ip->param->branch_strat == SOLVER_IP_SEARCH_DFS)
    {
        printf("\n");
        printf("Number of branch-and-bound nodes: %d\n",
               bac_env->ip->branch_count);
    }
    printf(" - Solution (Down) %.2f\n", cp->ip->lowerboundG);
    printf(" - Solution (Upper) %.2f\n", cp->ip->upperboundG);

    rval = 0;

CLEANUP:
    rval = stats_stop(bac_stats->total, !rval);
    printf("Total Running Time: %.2ld (milliseconds)",
           stats_get_current_time(bac_stats->total));
    printf("\n");

    return rval;
}
