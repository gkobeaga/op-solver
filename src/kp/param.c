#include "kp/kp.h"
#include "op-solver.h"

kp_param *
kp_create_param(void)
{
    kp_param *param      = malloc(sizeof(kp_param));
    param->time_limit    = 5 * 60 * 60 * 1000;
    param->exact_tech    = KP_EXACT_BAB;
    param->check_input   = 1;
    param->reorder_items = 0;
    param->epsilon       = 0.0001;
    param->p_epsilon     = 0.0001;
    param->w_epsilon     = 0.0001;
    return param;
}

void
kp_free_param(kp_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
