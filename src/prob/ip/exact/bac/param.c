#include "ip/exact/bac/bac.h"
#include "op-solver.h"

ip_exact_bac_param *
ip_create_exact_bac_param(void)
{
    ip_exact_bac_param *param = malloc(sizeof(ip_exact_bac_param));
    param->time_limit         = 5 * 60 * 60 * 1000;
    param->branch_strat       = SOLVER_IP_SEARCH_DFS;
    param->branch_select      = SOLVER_IP_SELECT_SIMPLE;
    // Should not be greater than this
    param->pruning_tol = 1.0 - SOLVER_PRICE_MAXPENALTY;
    return param;
}

void
ip_free_exact_bac_param(ip_exact_bac_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
