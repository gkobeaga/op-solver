#include "cp/cp.h"
#include "cp/init/init.h"
#include "cp/heur/heur.h"
#include "cp/exact/exact.h"

cp_env *
cp_create_env(void)
{
    cp_env *env    = malloc(sizeof(cp_env));
    env->verbosity = SOLVER_VERBOSITY_OFF;
    env->param     = cp_create_param();
    env->stats     = cp_create_stats();
    env->init      = cp_create_init_env();
    env->heur      = cp_create_heur_env();
#if HAVE_LP_SOLVER
    env->exact = cp_create_exact_env();
#else
    env->exact = NULL;

#endif
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
#if HAVE_LP_SOLVER
        cp_free_exact_env(&((*env)->exact));
#endif
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

#if HAVE_LP_SOLVER
    rval = cp_parse_exact_args(argc, argv, env->exact);
#endif

    return rval;
}
