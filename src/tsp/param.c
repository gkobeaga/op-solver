#include "op-solver.h"
#include "tsp/tsp.h"

tsp_param *
tsp_create_param(void)
{
    tsp_param *param  = malloc(sizeof(tsp_param));
    param->time_limit = 5 * 60 * 60 * 1000;
    param->exact      = SOLVER_TSP_EXACT_BAC;
    param->heur       = SOLVER_TSP_HEUR_EA;
    return param;
}

void
tsp_free_param(tsp_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
