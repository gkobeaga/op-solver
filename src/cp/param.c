#include "cp/cp.h"
#include "op-solver.h"

cp_param *
cp_create_param(void)
{
    cp_param *param   = malloc(sizeof(cp_param));
    param->time_limit = 5 * 60 * 60 * 1000;
#if HAVE_LP_SOLVER
    param->appr = SOLVER_CP_APPR_EXACT_BAC;
#else
    param->appr = SOLVER_CP_APPR_HEUR_EA;
#endif
    return param;
}

void
cp_free_param(cp_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
