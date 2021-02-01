#include "op-solver.h"
#include "cp/init/init.h"

cp_init_env *
cp_create_init_env(void)
{
    cp_init_env *env = malloc(sizeof(cp_init_env));
    env->verbosity   = SOLVER_VERBOSITY_OFF;
    env->param       = cp_create_init_param();
    env->stats       = cp_create_init_stats();
    return env;
}

void
cp_free_init_env(cp_init_env **env)
{
    if (*env)
    {
        cp_free_init_param(&(*env)->param);
        cp_free_init_stats(&(*env)->stats);
        free(*env);
        *env = NULL;
    }
}
