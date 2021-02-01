#include "cp/heur/ea/ea.h"
#include "op-solver.h"

cp_heur_ea_param *
cp_create_heur_ea_param(void)
{
    cp_heur_ea_param *param = malloc(sizeof(cp_heur_ea_param));
    param->time_limit       = 5 * 60 * 60 * 1000;
    param->it_lim           = INT_MAX;
    param->pop_size         = 100;
    param->pop_stop         = 25;
    param->d2d              = 50;
    param->nparsel          = 10;
    param->pmut             = 0.01;
    param->len_improve1     = 1;
    param->len_improve2     = 0;
    return param;
}

void
cp_free_heur_ea_param(cp_heur_ea_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
