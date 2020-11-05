#include "op-solver.h"

cp_exact_bac_env *
cp_create_exact_bac_env(void)
{
    cp_exact_bac_env *env          = malloc(sizeof(cp_exact_bac_env));
    env->verbosity                 = SOLVER_VERBOSITY_OFF;
    env->heur                      = cp_create_heur_env();
    env->heur->ea->param->pop_size = 10;
    env->heur->ea->param->d2d      = 5;
    env->heur->ea->param->nparsel  = 3;
    env->heur->ea->verbosity       = SOLVER_VERBOSITY_OFF;
    env->ip                        = ip_create_exact_bac_env();
    env->ip->orig_env              = env;
    env->ip->sep_loop              = cp_sep_loop;
    env->ip->xheur                 = cp_get_xheur_ea__;
    env->ip->pricing_loop          = cp_add_badvars;
    env->ip->recover_infeas        = cp_recover_infeas;
    env->ip->dual_bound            = cp_get_branch_dual_bound;
    env->ip->verify_infeas         = cp_verify_branch_infeas;
    env->ip->verify_prune          = cp_verify_branch_prune;
    env->ip->check_sol             = cp_is_ipsol_integral_connected;
    env->ip->find_branch           = cp_find_branch;
    env->kp                        = kp_create_env();
    env->stats                     = cp_create_exact_bac_stats();
    env->param                     = cp_create_exact_bac_param();
    env->cuts                      = NULL;
    env->cut_queue                 = NULL;
    return env;
}

void
cp_free_exact_bac_env(cp_exact_bac_env **env)
{
    if (*env)
    {
        while ((*env)->cut_queue)
        {
            cp_cut *cut       = (*env)->cut_queue;
            (*env)->cut_queue = cut->next;
            cp_free_cut(&cut);
        }
        cp_free_heur_env(&(*env)->heur);
        cp_free_exact_bac_stats(&(*env)->stats);
        cp_free_exact_bac_param(&(*env)->param);
        ip_free_exact_bac_env(&(*env)->ip);
        kp_free_env(&(*env)->kp);
        if ((*env)->cuts)
            cp_free_cut_repo(&(*env)->cuts);
        free(*env);
        *env = NULL;
    }
}

int
cp_parse_exact_bac_args(int argc, char *argv[], cp_exact_bac_env *env)
{
    int k;
    cp_exact_bac_param *param = env->param;

    for (k = 2; k < argc; k++)
    {

        if (!strcmp(argv[k], "--op-exact-bac-verbose"))
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
                env->verbosity = SOLVER_VERBOSITY_OFF;
                break;
            case 1:
                env->verbosity     = SOLVER_VERBOSITY_INFO;
                env->ip->verbosity = SOLVER_VERBOSITY_INFO;
                break;
            default:
                printf("Invalid argument value\n");
                return 1;
            }
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sep-sec-cc-2"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->sep_sec_cc_2 = atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sec-cc-srk"))
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
                param->srk_rule         = CP_SRK_NONE;
                param->srk_s2           = 0;
                param->srk_s3           = 0;
                param->srk_extra        = 0;
                param->sep_sec_cc_extra = 0;
                break;
            case 1:
                param->srk_rule  = CP_SRK_C1C2;
                param->srk_s2    = 0;
                param->srk_s3    = 1;
                param->srk_extra = 0;
                break;
            case 2:
                param->srk_rule  = CP_SRK_S1;
                param->srk_s2    = 0;
                param->srk_s3    = 1;
                param->srk_extra = 0;
                break;
            case 3:
                param->srk_rule  = CP_SRK_C1C2;
                param->srk_s2    = 0;
                param->srk_s3    = 1;
                param->srk_extra = 1;
                break;
            case 4:
                param->srk_rule  = CP_SRK_S1;
                param->srk_s2    = 0;
                param->srk_s3    = 1;
                param->srk_extra = 1;
                break;
            default:
                printf("Invalid argument value\n");
                return 1;
            }
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sep-sec-cc-extra"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->sep_sec_cc_extra = atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sep-blossom-fst"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->sep_blossom_mst = atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sep-blossom-eph"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->sep_blossom_fast = atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sep-blossom-egh"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->sep_blossom_ghfast = atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sep-cover-edge"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->sep_cover_edge = atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sep-cover-vertex"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->sep_cover_vertex = atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sep-cover-cycle"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->sep_cover_cycle = atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sep-path"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->sep_path = atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--op-exact-bac-sep-loop"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->sep_loop = atoi(argv[k]);
        }
        else if (!strcmp(argv[k], "--op-exact-bac-xheur-vph"))
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
                env->ip->xheur   = cp_get_xheur_greedy__;
                param->xheur_vph = 0;
                break;
            case 1:
                env->ip->xheur   = cp_get_xheur_ea__;
                param->xheur_vph = 1;
                break;
            default:
                printf("Invalid argument value\n");
                return 1;
            }
        }
        else if (!strcmp(argv[k], "--op-exact-bac-xheur-vph-meta"))
        {
            k++;
            if (k == argc || argv[k][0] == '\0' || argv[k][0] == '-')
            {
                printf("Missing argument value\n");
                return 1;
            }

            param->xheur_vph_meta = atoi(argv[k]);
        }
    }

    return 0;
}
