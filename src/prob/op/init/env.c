#include "op-solver.h"
#include "op/init/init.h"

op_init_env *
op_create_init_env(void)
{
    op_init_env *env = malloc(sizeof(op_init_env));
    env->verbosity   = SOLVER_VERBOSITY_OFF;
    return env;
}

void
op_free_init_env(op_init_env **env)
{
    if (*env)
    {
        free(*env);
        *env = NULL;
    }
}
