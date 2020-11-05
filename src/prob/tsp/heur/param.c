#include "op-solver.h"
#include "tsp/heur/heur.h"

tsp_heur_param *
tsp_create_heur_param(void)
{
    tsp_heur_param *param     = malloc(sizeof(tsp_heur_param));
    param->time_limit         = 5 * 60 * 60 * 1000;
    param->strat              = SOLVER_TSP_HEUR_LINKERN;
    param->two_and_a_half_opt = 0;
    return param;
}

void
tsp_free_heur_param(tsp_heur_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
