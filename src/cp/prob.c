#include "cp/cp.h"
#include "op-solver.h"

static int
eval_sol_obj(cp_prob *cp, cp_sol *sol);

static void
cp_init_prob_work(cp_prob *cp, solver_data *data)
{
    cp->n          = data->map->img_n;
    cp->cap        = data->cap;
    cp->data       = data;
    cp->sol_status = SOLVER_UNDEF;
    cp->sol        = cp_create_sol(cp);
#if HAVE_LP_SOLVER
    cp->ip = ip_create_prob();
#else
    cp->ip = NULL;
#endif
    cp->eval_sol_obj = eval_sol_obj;
    return;
}

cp_prob *
cp_create_prob(solver_data *data)
{
    cp_prob *cp = malloc(sizeof(cp_prob));
    cp_init_prob_work(cp, data);
    return cp;
}

static void
cp_delete_prob_work(cp_prob *cp);

void
cp_erase_prob(cp_prob *cp)
{
    solver_data *data = cp->data;
    cp_delete_prob_work(cp);
    cp_init_prob_work(cp, data);
    return;
}

static void
cp_delete_prob_work(cp_prob *cp)
{
#if HAVE_LP_SOLVER
    ip_free_prob(&(cp->ip));
#endif
    cp_free_sol(&(cp->sol));
    return;
}

void
cp_free_prob(cp_prob **cp)
{
    if (*cp)
    {
        cp_delete_prob_work(*cp);
        free(*cp);
        *cp = NULL;
    }
    return;
}

static int
eval_sol_obj(cp_prob *cp, cp_sol *sol)
{
    sol->val = 0.0;
    for (int i = 0; i < cp->n; i++)
    {
        if (sol->selected[i])
            sol->val += cp->data->obj_node[i];
    }
    return 0;
}
