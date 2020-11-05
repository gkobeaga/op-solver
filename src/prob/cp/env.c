#include "op-solver.h"

cp_env *
cp_create_env(void)
{
    cp_env *env    = malloc(sizeof(cp_env));
    env->verbosity = SOLVER_VERBOSITY_OFF;
    env->param     = cp_create_param();
    env->stats     = cp_create_stats();
    env->init      = cp_create_init_env();
    env->heur      = cp_create_heur_env();
    env->exact     = cp_create_exact_env();
    strcpy(env->sol_file, "prob.sol");
    return env;
}

void
cp_free_env(cp_env **env)
{
    if (*env)
    {
        cp_free_param(&((*env)->param));
        cp_free_stats(&((*env)->stats));
        cp_free_init_env(&((*env)->init));
        cp_free_heur_env(&((*env)->heur));
        cp_free_exact_env(&((*env)->exact));
        free(*env);
        *env = NULL;
    }
}

int
cp_parse_args(int argc, char *argv[], cp_env *env)
{
    int rval        = 0;
    cp_param *param = env->param;

    for (int k = 2; k < argc; k++)
    {

        if (!strcmp(argv[k], "--op-exact"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            switch (atoi(argv[k]))
            {
            case 0:
                param->appr = SOLVER_CP_APPR_HEUR_EA;
                break;
            case 1:
                param->appr = SOLVER_CP_APPR_EXACT_BAC;
                break;
            default:
                printf("Invalid argument value\n");
                return 1;
            }
        }
    }

    rval = cp_parse_exact_args(argc, argv, env->exact);
    return rval;
}
