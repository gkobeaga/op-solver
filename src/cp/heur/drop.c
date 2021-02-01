#include "cp/cp.h"
#include "op-solver.h"

int
cp_recover_heur_infeas(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol)
{
    int rval = 0;
    double cap, worst, cost, greedyeval;
    int i, j, tprev, tnext, firstit;
    int nodesel, nodeprev, nodenext;
    double *tcost;

    if (cp->n < 4)
    {
        fprintf(stderr, "Cannot drop nodes in an %d node tour\n", cp->n);
        return 1;
    }

    cap = sol->cap;

    tcost = malloc(cp->n * sizeof(double));

    firstit = 1;
    while (cap > cp->cap)
    {
        if (firstit)
        {
            for (i = 0; i < cp->n; i++)
            {
                if (sol->selected[i])
                {
                    tprev = sol->cod_bk[i];
                    tnext = sol->cod_fr[i];

                    tcost[i] = (double)(data_get_norm(cp->data, tprev, i) +
                                        data_get_norm(cp->data, i, tnext) -
                                        data_get_norm(cp->data, tprev, tnext));
                }
            }
            firstit = 0;
        }

        worst = SOLVER_MAXDOUBLE;
        for (i = 0; i < cp->n; i++)
        {
            if (sol->selected[i] && (i != cp->data->from || i != cp->data->to))
            {
                if (cp->data->obj_node[i] != 0)
                    greedyeval = cp->data->obj_node[i] / tcost[i];
                else
                    greedyeval = -SOLVER_MAXDOUBLE;

                if (greedyeval < worst)
                {
                    nodesel  = i;
                    nodeprev = sol->cod_bk[i];
                    nodenext = sol->cod_fr[i];
                    worst    = greedyeval;
                    cost     = tcost[i];
                }
            }
        }

        if (worst == SOLVER_MAXDOUBLE)
        {
            worst = 0.0;
            for (i = 0; i < cp->n; i++)
            {
                if (sol->selected[i] &&
                    (i != cp->data->from || i != cp->data->to))
                {
                    if (data_get_norm(cp->data, 0, i) > worst)
                    {
                        nodesel  = i;
                        nodeprev = sol->cod_bk[i];
                        nodenext = sol->cod_fr[i];
                        worst    = data_get_norm(cp->data, 0, i);
                        cost     = tcost[i];
                    }
                }
            }
        }

        sol->cod_fr[nodeprev] = nodenext;
        sol->cod_bk[nodenext] = nodeprev;
        sol->cod_fr[nodesel]  = nodesel;
        sol->cod_bk[nodesel]  = nodesel;

        tcost[nodeprev] =
        (double)(data_get_norm(cp->data, sol->cod_bk[nodeprev], nodeprev) +
                 data_get_norm(cp->data, nodeprev, nodenext) -
                 data_get_norm(cp->data, sol->cod_bk[nodeprev], nodenext));
        tcost[nodenext] =
        (double)(data_get_norm(cp->data, nodeprev, nodenext) +
                 data_get_norm(cp->data, nodenext, sol->cod_fr[nodenext]) -
                 data_get_norm(cp->data, nodeprev, sol->cod_fr[nodenext]));

        sol->selected[nodesel] = 0;
        sol->ns -= 1;
        sol->val -= cp->data->obj_node[nodesel];

        cap -= cost;
    }

    tprev = cp->data->from;
    cost  = 0;
    for (i = 0; i < sol->ns; i++)
    {
        sol->cycle[i] = tprev;
        cost += data_get_norm(cp->data, tprev, sol->cod_fr[tprev]);
        tprev = sol->cod_fr[tprev];
    }
    check_assert(cost == cap, "", CLEANUP);
    check_assert(tprev == cp->data->from, "", CLEANUP);
    for (i = sol->ns; i < cp->n; i++) sol->cycle[i] = -1;

    sol->cap = cap;
    free(tcost);

CLEANUP:

    return rval;
}

int
cp_drop_sol_node(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol, int node)
{
    int rval = 0;
    int i, j, prev, next;
    double cost;

    if (cp->n < 4)
    {
        fprintf(stderr, "Cannot drop nodes in an %d node tour\n", cp->n);
    }

    prev = sol->cod_bk[node];
    next = sol->cod_fr[node];

    cost = (double)(data_get_norm(cp->data, prev, node) +
                    data_get_norm(cp->data, node, next) -
                    data_get_norm(cp->data, prev, next));

    sol->cod_fr[prev] = next;
    sol->cod_bk[next] = prev;
    sol->cod_fr[node] = node;
    sol->cod_bk[node] = node;

    sol->selected[node] = 0;
    sol->ns -= 1;

    prev = 0;
    for (i = 0; i < sol->ns; i++)
    {
        sol->cycle[i] = prev;
        prev          = sol->cod_fr[prev];
    }
    assert(prev == 0);
    for (i = sol->ns; i < cp->n; i++) sol->cycle[i] = -1;

    j = 0;
    for (i = 0; i < cp->n; i++)
    {
        if (sol->selected[i])
        {
            sol->sposition[j] = i;
            j++;
        }
    }

    sol->cap -= cost;
    sol->val -= cp->data->obj_node[node];
    return rval;
}
