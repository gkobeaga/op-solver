#include "op-solver.h"
#include "op/op.h"

op_param *
op_create_param(void)
{
    op_param *param   = malloc(sizeof(op_param));
    param->time_limit = 5 * 60 * 60 * 1000;
    param->exact      = SOLVER_OP_EXACT_BAC;
    param->heur       = SOLVER_OP_HEUR_EA;
    return param;
}

void
op_free_param(op_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
