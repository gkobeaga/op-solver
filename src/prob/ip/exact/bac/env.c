#include "ip/exact/bac/bac.h"
#include "op-solver.h"

static int
undef_sep_loop(void *prob, void *env, int *nadded),
undef_pricing_loop(void *prob, void *env, int *nadded),
undef_verify_prune(void *prob, void *env, int *prune),
undef_verify_infeas(void *prob, void *env, int *prune),
undef_recover_infeas(void *prob, void *env),
undef_dual_bound(void *prob, void *env, mpf_t bound),
undef_xheur(void *prob, void *env, void *sol, int *nadded),
undef_check_integrality(ip_sol *sol),
undef_update_sol(void *prob, void *env, int *nadded),
undef_find_branch(void *prob, void *env, ip_branch **branch);

ip_exact_bac_env *
ip_create_exact_bac_env(void)
{
    ip_exact_bac_env *env = malloc(sizeof(ip_exact_bac_env));
    env->verbosity        = SOLVER_VERBOSITY_OFF;
    env->stats            = ip_create_exact_bac_stats();
    env->param            = ip_create_exact_bac_param();
    env->sep_loop         = undef_sep_loop;
    env->pricing_loop     = undef_pricing_loop;
    env->verify_prune     = undef_verify_prune;
    env->verify_infeas    = undef_verify_infeas;
    env->dual_bound       = undef_dual_bound;
    env->xheur            = undef_xheur;
    env->check_sol        = undef_check_integrality;
    env->update_sol       = undef_update_sol;
    env->find_branch      = undef_find_branch;

    env->ngot           = 0;
    env->depth          = 0;
    env->branch_log     = "branch.log";
    env->branch_count   = 0;
    env->history_space  = 50;
    env->history_depth  = 0;
    env->history        = malloc(env->history_space * sizeof(ip_branch *));
    env->lp             = lp_create_prob();
    env->branch_log     = NULL;
    env->id             = -1;
    env->parent_id      = -1;
    env->root           = 0;
    env->upperboundN    = SOLVER_MAXDOUBLE;
    env->lowerboundN    = -SOLVER_MAXDOUBLE;
    env->infeas_bound   = 0.0;
    env->farkas_pricing = 0;

    return env;
}

void
ip_free_exact_bac_env(ip_exact_bac_env **env)
{
    if (*env)
    {
        for (int i = 0; i < (*env)->history_depth; i++)
        {
            free((*env)->history[i]);
        }
        free((*env)->history);
        lp_free_prob(&(*env)->lp);
        ip_free_exact_bac_stats(&(*env)->stats);
        ip_free_exact_bac_param(&(*env)->param);
        free(*env);
        *env = NULL;
    }
}

static int
undef_sep_loop(void *prob, void *env, int *nadded)
{
    printf("Undefined sep_loop\n");
    return 1;
}

static int
undef_pricing_loop(void *prob, void *env, int *nadded)
{
    printf("Undefined pricing_loop\n");
    return 1;
}
static int
undef_verify_prune(void *prob, void *env, int *prune)
{
    printf("Undefined verify_prune\n");
    return 1;
}
static int
undef_verify_infeas(void *prob, void *env, int *prune)
{
    printf("Undefined verify_infeas\n");
    return 1;
}
static int
undef_recover_infeas(void *prob, void *env)
{
    printf("Undefined recover_infeas\n");
    return 1;
}
static int
undef_dual_bound(void *prob, void *env, mpf_t bound)
{
    printf("Undefined dual_bound\n");
    return 1;
}
static int
undef_xheur(void *prob, void *env, void *sol, int *nadded)
{
    printf("Undefined xheur\n");
    return 1;
}
static int
undef_check_integrality(ip_sol *sol)
{
    printf("Undefined check_integrality\n");
    return 0;
}
static int
undef_update_sol(void *prob, void *env, int *nadded)
{
    printf("Undefined update_sol\n");
    return 1;
}
static int
undef_find_branch(void *prob, void *env, ip_branch **branch)
{
    printf("Undefined find_branch\n");
    return 1;
}
