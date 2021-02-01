#include "cp/heur/heur.h"
#include "op-solver.h"

cp_heur_param *
cp_create_heur_param(void)
{
    cp_heur_param *param = malloc(sizeof(cp_heur_param));
    param->time_limit    = 5 * 60 * 60 * 1000;
    param->improve_sol   = (SOLVER_CP_ADD_SD | SOLVER_CP_ADD_3N);
    param->pinit         = 0.5;
    param->init          = SOLVER_CP_INIT_RAND;
    param->select        = SOLVER_CP_SEL_BERNOULLI;
    return param;
}

void
cp_free_heur_param(cp_heur_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
