#include "op-solver.h"
#include "op/init/init.h"

op_init_param *
op_create_init_param(void)
{
    op_init_param *param = malloc(sizeof(op_init_param));
    param->time_limit    = 5 * 60 * 60 * 1000;
    return param;
}

void
op_free_init_param(op_init_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
