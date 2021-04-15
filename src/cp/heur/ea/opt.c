#include "cp/heur/ea/ea.h"
#include "cp/cp.h"
#include "op-solver.h"

int
cp_opt_heur_ea(cp_prob *cp, cp_heur_env *heur_env, cp_pop *pop, cp_sol *sol)
{
    int rval = 0;
    int it, i;
    cp_heur_ea_env *ea_env     = heur_env->ea;
    cp_heur_ea_param *ea_param = ea_env->param;
    cp_heur_ea_stats *ea_stats = ea_env->stats;
    int *parent                = NULL;
    cp_sol *child              = NULL;
    cp_sol *best_sol;
    cp_sol *pop_sol;
    int stop_ind;
    double stop_val;
    int stop         = 0;
    int sol_improved = 0;

    rval = stats_start(ea_stats->total);
    check_rval(rval, "stats_start failed", cleanup);

    parent = malloc(2 * sizeof(int));
    child  = cp_create_sol(cp);

    if (cp->n < 4)
    {
        rval = stats_stop(ea_stats->total, 0);
        check_rval(rval, "stats_stop failed", cleanup);
        printf("Less than 4 node problem\n");
        goto done;
    }

    for (it = 1; (it < ea_param->it_lim + 1) && !stop; it++)
    {
        sol_improved = 0;
        rval         = stats_start(ea_stats->it);
        check_rval(rval, "stats_start failed", cleanup);
        if (it % ea_param->d2d != 0)
        {
            rval = cp_heur_ea_selection(cp, ea_env, pop, parent);
            check_rval(rval, "failed", cleanup);

            if (parent[0] != parent[1])
            {
                cp_heur_ea_crossover(cp, ea_env, pop->sol[parent[0]],
                                     pop->sol[parent[1]], child);
            }
            else
            {
                cp_copy_sol(pop->sol[parent[0]], child);
            }

            if (rng_bernoulli(ea_param->pmut))
                cp_heur_ea_mutation(cp, ea_env, child);

            if (pop->best_val < child->val)
            {
                sol_improved = 1;
            }
            if (pop->worst_val < child->val)
            {
                cp_set_pop_sol(pop, child, pop->worst_ind);
                cp_update_pop(pop);
            }
            cp_erase_sol(child);
        }
        else
        {
            rval = stats_start(ea_stats->infeas_recover);
            check_rval(rval, "stats_start failed", cleanup);
            if (ea_param->len_improve1)
            {
                for (i = 0; i < pop->size; i++)
                {
                    pop_sol = pop->sol[i];
                    cp->eval_sol_obj(cp, pop_sol);
                    cp_improve_heur_cycle_length(cp, heur_env, pop_sol);
                }
            }

            if (cp->cap != 0.0)
            {
                for (i = 0; i < pop->size; i++)
                {
                    pop_sol = pop->sol[i];
                    heur_env->recover_infeas(cp, heur_env, pop_sol);
                    heur_env->local_search(cp, heur_env, pop_sol);
                }
            }

            if (ea_param->len_improve2)
            {
                for (i = 0; i < pop->size; i++)
                {
                    pop_sol = pop->sol[i];
                    cp_improve_heur_cycle_length(cp, heur_env, pop_sol);
                }
            }

            sol_improved = 1;

            cp_update_pop(pop);
            best_sol = pop->sol[pop->best_ind];
            if (sol->val < best_sol->val)
                cp_copy_sol(best_sol, sol);

            rval = stats_stop(ea_stats->infeas_recover, 1);
            check_rval(rval, "stats_stop failed", cleanup);
        }

        if (it % ea_param->d2d == 0 || cp->cap == 0.0)
        {
            if (ea_env->verbosity >= SOLVER_VERBOSITY_INFO)
                printf("cp   | EA :  %d it : best %.0f : worst %.0f (%.2ld "
                       "millisec) \n",
                       it, pop->best_val, pop->worst_val,
                       stats_get_current_time(ea_stats->total));

            if (ea_param->pop_stop)
            {

                stop_ind = floor(pop->size / 100.0 * ea_param->pop_stop);
                stop_val = pop->sol[stop_ind]->val;
                if (pop->best_val == stop_val)
                    stop = 1;
            }

            if (stats_get_current_time(ea_stats->total) > ea_param->time_limit)
            {
                stop = 1;
            }
        }

        rval = stats_stop(ea_stats->it, sol_improved);
        check_rval(rval, "stats_stop failed", cleanup);
    }

    free(parent);
    cp_free_sol(&child);

done:

    if (ea_stats->write_stats)
        cp_write_heur_ea_stats(cp, ea_env);

cleanup:

    rval = stats_stop(ea_stats->total, 1);
    check_rval(rval, "stats_stop failed", cleanup);
    if (sol->val > cp->sol->val)
        cp_copy_sol(sol, cp->sol);
    cp->sol_status = SOLVER_FEAS;

    return rval;
}
