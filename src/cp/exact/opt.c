#include "cp/cp.h"
#include "cp/exact/bac/bac.h"

int
cp_opt_exact(cp_prob *cp, cp_exact_env *exact_env, cp_sol *sol)
{
    int rval                    = 0;
    cp_exact_stats *exact_stats = exact_env->stats;
    cp_exact_param *exact_param = exact_env->param;

    rval = stats_start(exact_stats->total);
    check_rval(rval, "stats_start failed", CLEANUP);

    if (exact_param->appr == SOLVER_CP_APPR_EXACT_BAC)
    {
        rval = cp_opt_exact_bac(cp, exact_env->bac, sol);
        check_rval(rval, "failed", CLEANUP);
    }

    rval = 0;

CLEANUP:
    stats_stop(exact_stats->total, !rval);

    return rval;
}
