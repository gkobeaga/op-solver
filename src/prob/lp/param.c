#include "op-solver.h"
#include "prob/lp/lp.h"

lp_param *
lp_create_param(void)
{
    lp_param *param   = malloc(sizeof(lp_param));
    param->init_phase = 1;
    return param;
}

void
lp_free_param(lp_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
