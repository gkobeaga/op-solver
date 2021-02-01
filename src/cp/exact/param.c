#include "cp/exact/exact.h"
#include "op-solver.h"

cp_exact_param *
cp_create_exact_param(void)
{
    cp_exact_param *param = malloc(sizeof(cp_exact_param));
    param->time_limit     = 5 * 60 * 60 * 1000;
    param->appr           = SOLVER_CP_APPR_EXACT_BAC;
    return param;
}

void
cp_free_exact_param(cp_exact_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
