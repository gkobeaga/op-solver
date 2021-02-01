#include "ip/exact/exact.h"
#include "op-solver.h"

ip_exact_env *
ip_create_exact_env(void)
{
    ip_exact_env *env = malloc(sizeof(ip_exact_env));
    env->verbosity    = SOLVER_VERBOSITY_OFF;
    return env;
}

void
ip_free_exact_env(ip_exact_env **env)
{
    if (*env)
    {
        free(*env);
        *env = NULL;
    }
}
