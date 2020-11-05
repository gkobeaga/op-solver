#include "op-solver.h"
#include "cp/cp.h"
#include "cp/init/init.h"

static void
select_bernoulli(cp_prob *cp, cp_heur_env *env, cp_sol *sol);

int
cp_select_sol_vertices(cp_prob *cp, cp_heur_env *env, cp_sol *sol)
{
    int rval             = 0;
    cp_heur_param *param = env->param;

    if (param->select == SOLVER_CP_SEL_BERNOULLI)
        select_bernoulli(cp, env, sol);

    return rval;
}

static void
select_bernoulli(cp_prob *cp, cp_heur_env *env, cp_sol *sol)
{
    int i, ns;
    int from             = cp->data->from;
    int to               = cp->data->to;
    cp_heur_param *param = env->param;

    ns = 0;
    do
    {
        for (i = 0; i < cp->n; i++)
        {
            if (!sol->selected[from] && i == from)
            {
                sol->selected[from] = 1;
                sol->cycle[ns++]    = from;
            }
            else if (!sol->selected[to] && i == to)
            {
                sol->selected[to] = 1;
                sol->cycle[ns++]  = to;
            }
            else if (!sol->selected[i])
            {
                if (drand48() < param->pinit)
                {
                    sol->selected[i] = 1;
                    sol->cycle[ns++] = i;
                }
                else if (!sol->selected[i])
                    sol->selected[i] = 0;
            }
        }
    } while (ns <= 3);
    sol->ns = ns;

    int prev, first;
    sol->val = 0.0;
    sol->cap = 0.0;
    prev     = -1;
    for (i = 0; i < cp->n; i++)
    {
        if (sol->selected[i])
        {
            if (prev >= 0)
            {
                sol->cod_fr[prev] = i;
                sol->cod_bk[i]    = prev;
                sol->cap += data_get_norm(cp->data, prev, i);
            }
            else
                first = i;
            prev = i;
            sol->val += cp->data->obj_node[i];
        }
        else
        {
            sol->cod_fr[i] = i;
            sol->cod_bk[i] = i;
        }
    }
    sol->cod_bk[first] = prev;
    sol->cod_fr[prev]  = first;
    sol->cap += data_get_norm(cp->data, prev, first);
    return;
}

int
cp_select_best3nodes(cp_prob *cp, cp_heur_env *env, cp_sol *sol)
{
    int rval = 0;
    int i, j, v1, v2;
    double best_len, best_val, val_tmp, cap_tmp;

    best_len = -SOLVER_MAXDOUBLE;
    best_val = -SOLVER_MAXDOUBLE;
    for (i = 1; i < cp->n; i++)
    {
        for (j = i + 1; j < cp->n; j++)
        {
            val_tmp = cp->data->obj_node[0] + cp->data->obj_node[i] +
                      cp->data->obj_node[j];
            cap_tmp = data_get_norm(cp->data, 0, i) +
                      data_get_norm(cp->data, 0, j) +
                      data_get_norm(cp->data, i, j);
            if (cap_tmp < cp->cap && val_tmp > best_val)
            {
                best_val = val_tmp;
                best_len = cap_tmp;
                v1       = i;
                v2       = j;
            }
            else if (cap_tmp < best_len && val_tmp == best_val)
            {
                best_val = val_tmp;
                best_len = cap_tmp;
                v1       = i;
                v2       = j;
            }
        }
    }

    if (best_val == -SOLVER_MAXDOUBLE)
    {
        rval = 1;
        goto DONE;
    }

    sol->selected[0]  = 1;
    sol->selected[v1] = 1;
    sol->selected[v2] = 1;

    for (i = 0; i < cp->n; i++)
    {
        sol->cod_fr[i] = i;
        sol->cod_bk[i] = i;
    }

    sol->cod_fr[0]  = v1;
    sol->cod_bk[v1] = 0;
    sol->cod_fr[v1] = v2;
    sol->cod_bk[v2] = v1;
    sol->cod_fr[v2] = 0;
    sol->cod_bk[0]  = v2;

    sol->cycle[0] = 0;
    sol->cycle[1] = v1;
    sol->cycle[2] = v2;
    sol->ns       = 3;
    sol->val      = best_val;
    sol->cap      = best_len;

DONE:
    return rval;
}
