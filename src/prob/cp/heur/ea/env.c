#include "cp/heur/ea/ea.h"
#include "op-solver.h"

cp_heur_ea_env *
cp_create_heur_ea_env(void)
{
    cp_heur_ea_env *env = malloc(sizeof(cp_heur_ea_env));
    env->verbosity      = SOLVER_VERBOSITY_INFO;
    env->param          = cp_create_heur_ea_param();
    env->stats          = cp_create_heur_ea_stats();
    return env;
}

void
cp_free_heur_ea_env(cp_heur_ea_env **env)
{
    if (*env)
    {
        cp_free_heur_ea_param(&((*env)->param));
        cp_free_heur_ea_stats(&((*env)->stats));
        free(*env);
        *env = NULL;
    }
}
