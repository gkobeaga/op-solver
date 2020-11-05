#include "op-solver.h"
#include "tsp/tsp.h"
#include "tsp/init/init.h"

static void
randcycle(int ncount, int *cyc);

int
tsp_init_sol_random(tsp_prob *tsp, tsp_init_env *env, tsp_sol *sol)
{
    int i;
    if (env->verbosity >= SOLVER_VERBOSITY_INFO)
        printf("tsp   :  Generating Random Rour...");
    randcycle(tsp->n, sol->cycle);
    sol->val = data_get_norm(tsp->data, sol->cycle[tsp->n - 1], sol->cycle[0]);
    for (i = 1; i < tsp->n; i++)
    {
        sol->cod_fr[sol->cycle[i - 1]] = sol->cycle[i];
        sol->cod_bk[sol->cycle[i]]     = sol->cycle[i - 1];
        sol->val += data_get_norm(tsp->data, sol->cycle[i - 1], sol->cycle[i]);
    }
    sol->ns                             = tsp->n;
    sol->tot_n                          = tsp->n;
    sol->cod_fr[sol->cycle[tsp->n - 1]] = sol->cycle[0];
    sol->cod_bk[sol->cycle[0]]          = sol->cycle[tsp->n - 1];
    memset(sol->selected, 1, tsp->n * sizeof(int));
    if (env->verbosity >= SOLVER_VERBOSITY_INFO)
        printf("tsp   :  Generating Random Rour...");
    return 0;
}

static void
randcycle(int ncount, int *cyc)
{
    int i, k, temp;
    for (i = 0; i < ncount; i++) cyc[i] = i;

    for (i = ncount; i > 2; i--)
    {
        k = rand() % i;
        if (k != 0)
            SWAP(cyc[i - 1], cyc[k], temp);
    }
    return;
}
