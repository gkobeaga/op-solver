#include "op-solver.h"
#include "tsp/heur/linkern/linkern.h"

tsp_heur_linkern_env *
tsp_create_heur_linkern_env(void)
{
    tsp_heur_linkern_env *env = malloc(sizeof(tsp_heur_linkern_env));
    env->verbosity            = SOLVER_VERBOSITY_OFF;
    env->param                = tsp_create_heur_linkern_param();
    env->stats                = tsp_create_heur_linkern_stats();
    return env;
}

void
tsp_free_heur_linkern_env(tsp_heur_linkern_env **env)
{
    if (*env)
    {
        tsp_free_heur_linkern_param(&(*env)->param);
        tsp_free_heur_linkern_stats(&(*env)->stats);
        free(*env);
        *env = NULL;
    }
}
