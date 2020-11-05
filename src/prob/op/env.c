#include "op-solver.h"
#include "op/op.h"

op_env *
op_create_env(void)
{
    op_env *env    = malloc(sizeof(op_env));
    env->verbosity = SOLVER_VERBOSITY_OFF;
    env->param     = op_create_param();
    env->stats     = op_create_stats();
    env->init      = op_create_init_env();
    env->heur      = op_create_heur_env();
    env->exact     = op_create_exact_env();
    return env;
}

void
op_free_env(op_env **env)
{
    if (*env)
    {
        op_free_param(&((*env)->param));
        op_free_stats(&((*env)->stats));
        op_free_init_env(&((*env)->init));
        op_free_heur_env(&((*env)->heur));
        op_free_exact_env(&((*env)->exact));
        free(*env);
        *env = NULL;
    }
}
