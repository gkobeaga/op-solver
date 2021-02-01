#include "cp/heur/heur.h"
#include "tsp/heur/heur.h"

cp_heur_env *
cp_create_heur_env(void)
{
    cp_heur_env *env          = malloc(sizeof(cp_heur_env));
    env->verbosity            = SOLVER_VERBOSITY_OFF;
    env->param                = cp_create_heur_param();
    env->stats                = cp_create_heur_stats();
    env->ea                   = cp_create_heur_ea_env();
    env->tsp                  = tsp_create_heur_env();
    env->select_best3vertices = cp_select_best3nodes;
    env->recover_infeas       = cp_recover_heur_infeas;
    env->local_search         = cp_improve_heur_3n;
    return env;
}

void
cp_free_heur_env(cp_heur_env **env)
{
    if (*env)
    {
        cp_free_heur_param(&((*env)->param));
        cp_free_heur_stats(&((*env)->stats));
        cp_free_heur_ea_env(&(*env)->ea);
        tsp_free_heur_env(&(*env)->tsp);
        free(*env);
        *env = NULL;
    }
}

void
cp_conf_heur_env(cp_prob *cp, cp_heur_env *heur_env)
{
    cp_heur_param *heur_param = heur_env->param;

    heur_env->recover_infeas = cp_recover_heur_infeas;

    if (heur_param->improve_sol == SOLVER_CP_ADD_3N)
        heur_env->local_search = cp_improve_heur_3n;
    else if (heur_param->improve_sol == SOLVER_CP_ADD_IN)
        heur_env->local_search = cp_improve_heur_in;
}
