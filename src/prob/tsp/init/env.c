#include "op-solver.h"
#include "tsp/init/init.h"

tsp_init_env *
tsp_create_init_env(void)
{
    tsp_init_env *env = malloc(sizeof(tsp_init_env));
    env->verbosity    = SOLVER_VERBOSITY_OFF;
    env->param        = tsp_create_init_param();
    env->stats        = tsp_create_init_stats();
    return env;
}

void
tsp_free_init_env(tsp_init_env **env)
{
    if (*env)
    {
        tsp_free_init_param(&(*env)->param);
        tsp_free_init_stats(&(*env)->stats);
        free(*env);
        *env = NULL;
    }
}
