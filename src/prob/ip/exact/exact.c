#include "op-solver.h"
#include "ip/ip.h"

int
ip_opt_exact(ip_prob *ip, ip_exact_env *env, ip_sol *sol)
{
    int rval                    = 0;
    ip_exact_stats *exact_stats = env->stats;
    ip_exact_param *exact_param = env->param;

    rval = stats_start(exact_stats->total);
    check_rval(rval, "stats_start failed", CLEANUP);

    if (exact_param->appr == SOLVER_IP_EXACT_APPR_BAC)
    {
        ip_opt_exact_bac(ip, env->bac, sol);
    }

CLEANUP:
    rval = stats_stop(exact_stats->total, !rval);

    return rval;
}
