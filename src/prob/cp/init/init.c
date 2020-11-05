#include "cp/cp.h"
#include "op-solver.h"
#include "tsp/tsp.h"

int
cp_init_sol(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol)
{
    int rval             = 0;
    cp_heur_param *param = heur_env->param;

    if (param->init == SOLVER_CP_INIT_BEST3)
    {
        heur_env->select_best3vertices(cp, heur_env, sol);
    }
    else if (param->init == SOLVER_CP_INIT_RAND)
    {
        if (!sol->ns)
            cp_select_sol_vertices(cp, heur_env, sol);
        cp_improve_heur_cycle_length(cp, heur_env, sol);
        heur_env->recover_infeas(cp, heur_env, sol);
        cp_improve_heur_cycle_length(cp, heur_env, sol);
    }
    else
        rval = 1;

    heur_env->local_search(cp, heur_env, sol);

    return rval;
}

int
cp_init_pop(cp_prob *cp, cp_heur_env *env, cp_pop *pop)
{
    int rval = 0;
    int i;
    cp_heur_stats *heur_stats = env->stats;
    cp_heur_param *heur_param = env->param;

    if (env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("cp   : > Population size: %d\n", pop->size);
    }
    for (i = 0; i < pop->size; i++)
    {
        cp_sol *sol = pop->sol[i];
        cp_erase_sol(sol);
        cp_init_sol(cp, env, sol);
        if (env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf(" %d: nvis: %d, capacity %.0f, value %.0f\n", i, sol->ns,
                   sol->cap, sol->val);
        }
        if (0 &&
            stats_get_current_time(heur_stats->total) > heur_param->time_limit)
        {
            pop->size = env->param->pop_size = i + 1;
            break;
        }
    }
    cp_update_pop(pop);

    return rval;
}
