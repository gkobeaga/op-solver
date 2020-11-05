#include "op-solver.h"

int
cp_opt(cp_prob *cp, cp_env *env, cp_sol *sol)
{
    int rval = 0;

    cp_stats *stats = env->stats;
    cp_param *param = env->param;

    rval = stats_start(stats->total);
    check_rval(rval, "stats_start failed", CLEANUP);
    if (param->appr == SOLVER_CP_APPR_HEUR_EA)
    {
        cp_pop *pop = cp_create_pop(cp, env->heur->ea->param->pop_size);
        rval        = cp_init_pop(cp, env->heur, pop);
        rval        = cp_opt_heur_ea(cp, env->heur, pop, sol);
        cp_free_pop(&pop);
        check_rval(rval, "failed", CLEANUP);
    }
    else if (param->appr == SOLVER_CP_APPR_EXACT_BAC)
    {
        rval = cp_opt_exact(cp, env->exact, sol);
        check_rval(rval, "failed", CLEANUP);
    }

    rval = cp_write_sol(cp, sol, env->sol_file);
    check_rval(rval, "failed", CLEANUP);

    rval = 0;

CLEANUP:
    rval = stats_stop(stats->total, !rval);
    printf("Total Running Time: %.2ld (milliseconds)",
           stats_get_current_time(stats->total));
    printf("\n");

    return rval;
}
