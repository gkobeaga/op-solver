#include "kp/kp.h"
#include "op-solver.h"

kp_env *
kp_create_env(void)
{
    kp_env *env        = malloc(sizeof(kp_env));
    env->verbosity     = SOLVER_VERBOSITY_OFF;
    env->param         = kp_create_param();
    env->stats         = kp_create_stats();
    env->branch_cap    = 0.0;
    env->branch_profit = 0.0;
    env->j             = 0;
    env->r             = 0;
    env->p_integer     = 0;
    env->w_integer     = 0;
    env->counter       = 0;
    env->max_counter   = 100000;
    return env;
}

void
kp_free_env(kp_env **env)
{
    if (*env)
    {
        kp_free_param(&((*env)->param));
        kp_free_stats(&((*env)->stats));
        free(*env);
        *env = NULL;
    }
}
