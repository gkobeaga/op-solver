#include "cp/heur/ea/ea.h"
#include "op-solver.h"

int
cp_heur_ea_selection(cp_prob *cp, cp_heur_ea_env *ea_env, cp_pop *pop,
                     int *parents)
{
    int rval                   = 0;
    cp_heur_ea_param *ea_param = ea_env->param;
    int i;
    double min                 = SOLVER_MAXDOUBLE;
    double sum                 = 0;
    double *reproductive_probs = NULL, *tprobs = NULL;
    int *preselected = NULL, *selected = NULL;
    int *indexes = NULL;
    cp_sol *sol;

    reproductive_probs = malloc(ea_param->nparsel * sizeof(double));
    check_null(reproductive_probs, "out of memory", cleanup);
    tprobs = malloc(ea_param->nparsel * sizeof(double));
    check_null(tprobs, "out of memory", cleanup);
    selected = malloc(2 * sizeof(int));
    check_null(selected, "out of memory", cleanup);
    preselected = malloc(ea_param->nparsel * sizeof(int));
    check_null(preselected, "out of memory", cleanup);
    indexes = malloc(ea_param->pop_size * sizeof(int));
    check_null(indexes, "out of memory", cleanup);
    for (i = 0; i < ea_param->pop_size; i++) indexes[i] = i;
    rng_choose(indexes, ea_param->pop_size, preselected, ea_param->nparsel,
               sizeof(int));
    for (i = 0; i < ea_param->nparsel; i++)
    {
        sol = pop->sol[preselected[i]];
        if (sol->val < min)
            min = sol->val;
    }
    // We add 1 to ensure that are not null.
    for (i = 0; i < ea_param->nparsel; i++)
        tprobs[i] = pop->sol[preselected[i]]->val - min + 1;
    for (i = 0; i < ea_param->nparsel; i++) sum += tprobs[i];
    for (i = 0; i < ea_param->nparsel; i++)
        reproductive_probs[i] = tprobs[i] / sum;
    rng_choose_p(preselected, reproductive_probs, ea_param->nparsel, selected,
                 2, sizeof(int));

    parents[0] = selected[0];
    parents[1] = selected[1];
cleanup:
    free(reproductive_probs);
    free(tprobs);
    free(selected);
    free(preselected);
    free(indexes);
    return rval;
}
