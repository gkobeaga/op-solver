#include "ip/ip.h"
#include "op-solver.h"

ip_env *
ip_create_env(void)
{
    ip_env *env    = malloc(sizeof(ip_env));
    env->verbosity = SOLVER_VERBOSITY_OFF;
    env->param     = ip_create_param();
    env->stats     = ip_create_stats();
    return env;
}

void
ip_free_env(ip_env **env)
{
    if (*env)
    {
        ip_free_param(&((*env)->param));
        ip_free_stats(&((*env)->stats));
        free(*env);
        *env = NULL;
    }
}
