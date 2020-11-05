#include "op-solver.h"
#include "tsp/heur/linkern/linkern.h"

tsp_heur_linkern_param *
tsp_create_heur_linkern_param(void)
{
    tsp_heur_linkern_param *param = malloc(sizeof(tsp_heur_linkern_param));
    param->time_limit             = 5 * 60 * 60 * 1000;
    param->length_bound           = 0.0;
    param->nruns                  = 1;
    param->kick_type              = TSP_LINKERN_WALK_KICK;
    param->nkicks                 = 1;
    return param;
}

void
tsp_free_heur_linkern_param(tsp_heur_linkern_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
