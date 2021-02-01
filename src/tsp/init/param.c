#include "op-solver.h"
#include "tsp/init/init.h"

tsp_init_param *
tsp_create_init_param(void)
{
    tsp_init_param *param = malloc(sizeof(tsp_init_param));
    param->time_limit     = 5 * 60 * 60 * 1000;
    param->init           = TSP_INIT_RANDOM;
    return param;
}

void
tsp_free_init_param(tsp_init_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
