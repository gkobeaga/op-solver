#include "tsp/tsp.h"
#include "tsp/init/init.h"
#include "tsp/heur/heur.h"

tsp_env *
tsp_create_env(void)
{
    tsp_env *env   = malloc(sizeof(tsp_env));
    env->verbosity = SOLVER_VERBOSITY_OFF;
    env->param     = tsp_create_param();
    env->stats     = tsp_create_stats();
    env->init      = tsp_create_init_env();
    env->heur      = tsp_create_heur_env();
    return env;
}

void
tsp_free_env(tsp_env **env)
{
    if (*env)
    {
        tsp_free_param(&((*env)->param));
        tsp_free_stats(&((*env)->stats));
        tsp_free_init_env(&((*env)->init));
        tsp_free_heur_env(&((*env)->heur));
        free(*env);
        *env = NULL;
    }
}
