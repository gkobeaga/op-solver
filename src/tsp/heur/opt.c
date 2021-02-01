#include "tsp/heur/heur.h"
#include "op-solver.h"

int
tsp_opt_heur(tsp_prob *tsp, tsp_heur_env *env, tsp_sol *sol)
{
    int rval              = 0;
    tsp_heur_stats *stats = env->stats;
    tsp_heur_param *param = env->param;

    rval = stats_start(stats->total);
    check_rval(rval, "stats_start failed", done);

    if (tsp->n <= 3)
    {
        fprintf(stderr,
                "tsp    : Cannot run local search in an %d node graph\n",
                tsp->n);
        rval = 1;
        goto done;
    }

    switch (param->strat)
    {
    case SOLVER_TSP_HEUR_NO:
        if (env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("tsp   :  warning no local search.");
        goto done;
    case SOLVER_TSP_HEUR_TWOOPT:
        if (env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("2-opt not implemented\n");
        goto done;
    case SOLVER_TSP_HEUR_THREEOPT:
        if (env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("3-opt not implemented\n");
        goto done;
    case SOLVER_TSP_HEUR_LINKERN:

        if (tsp_opt_heur_linkern(tsp, env->linkern, sol))
        {
            fprintf(stderr, "tsp   :   call_threeopt_tour failed\n");
            rval = 1;
            goto done;
        }
        goto done;
    default:
        fprintf(stderr, "tsp   :   invalid local search flag.\n");
        rval = 1;
        goto done;
    }
done:

    rval = stats_stop(stats->total, rval);

    return rval;
}
