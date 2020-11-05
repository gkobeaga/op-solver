#include "cp/exact/exact.h"
#include "op-solver.h"

cp_exact_env *
cp_create_exact_env(void)
{
    cp_exact_env *env = malloc(sizeof(cp_exact_env));
    env->verbosity    = SOLVER_VERBOSITY_OFF;
    env->stats        = cp_create_exact_stats();
    env->param        = cp_create_exact_param();
    env->bac          = cp_create_exact_bac_env();
    return env;
}

void
cp_free_exact_env(cp_exact_env **env)
{
    if (*env)
    {
        cp_free_exact_stats(&(*env)->stats);
        cp_free_exact_param(&(*env)->param);
        cp_free_exact_bac_env(&(*env)->bac);
        free(*env);
        *env = NULL;
    }
}

int
cp_parse_exact_args(int argc, char *argv[], cp_exact_env *env)
{
    int rval = 0;
    rval     = cp_parse_exact_bac_args(argc, argv, env->bac);
    return rval;
}
