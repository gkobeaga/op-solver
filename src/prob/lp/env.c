#include "op-solver.h"
#include "prob/lp/lp.h"

lp_env *
lp_create_env(void)
{
    lp_env *env    = malloc(sizeof(lp_env));
    env->verbosity = SOLVER_VERBOSITY_OFF;
    env->param     = lp_create_param();
    env->stats     = lp_create_stats();
    return env;
}

void
lp_free_env(lp_env **env)
{
    if (*env)
    {
        lp_free_param(&((*env)->param));
        lp_free_stats(&((*env)->stats));
        free(*env);
        *env = NULL;
    }
}
