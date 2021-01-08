#include "cp/init/init.h"
#include "op-solver.h"

cp_init_param *
cp_create_init_param(void)
{
    cp_init_param *param = malloc(sizeof(cp_init_param));
    param->time_limit    = 5 * 60 * 60 * 1000;
    param->pinit         = 0;
    param->init          = SOLVER_CP_INIT_RAND;
    param->select        = SOLVER_CP_SEL_BERNOULLI;
    return param;
}

void
cp_free_init_param(cp_init_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
