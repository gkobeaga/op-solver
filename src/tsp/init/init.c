#include "op-solver.h"
#include "tsp/tsp.h"
#include "tsp/init/init.h"

int
tsp_init_sol(tsp_prob *tsp, tsp_init_env *env, tsp_sol *sol)
{
    int rval              = 0;
    tsp_init_stats *stats = env->stats;
    tsp_init_param *param = env->param;

    rval = stats_start(stats->total);
    check_rval(rval, "stats_start failed", done);

    switch (param->init)
    {
    case TSP_INIT_RANDOM:
        if (tsp_init_sol_random(tsp, env, sol))
        {
            fprintf(stderr, "tsp   :   tsp_init_sol_randcycle failed\n");
            tsp_free_sol(&sol);
            rval = 1;
        }
        goto done;
    default:
        fprintf(stderr, "tsp   :   invalid tsp sol initialization\n");
        tsp_free_sol(&sol);
        rval = 1;
    }

done:

    rval = stats_stop(stats->total, rval);

    return rval;
}
