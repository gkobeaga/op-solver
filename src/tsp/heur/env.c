#include "op-solver.h"
#include "tsp/heur/heur.h"

tsp_heur_env *
tsp_create_heur_env(void)
{
    tsp_heur_env *env = malloc(sizeof(tsp_heur_env));
    env->verbosity    = SOLVER_VERBOSITY_INFO;
    env->stats        = tsp_create_heur_stats();
    env->param        = tsp_create_heur_param();
    env->linkern      = tsp_create_heur_linkern_env();
    return env;
}

void
tsp_free_heur_env(tsp_heur_env **env)
{
    if (*env)
    {
        tsp_free_heur_param(&(*env)->param);
        tsp_free_heur_stats(&(*env)->stats);
        tsp_free_heur_linkern_env(&(*env)->linkern);
        free(*env);
        *env = NULL;
    }
}
