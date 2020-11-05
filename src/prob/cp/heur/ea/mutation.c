#include "op-solver.h"
#include "cp/heur/ea/ea.h"

int
cp_heur_ea_mutation(cp_prob *cp, cp_heur_ea_env *ea_env, cp_sol *sol)
{
    int rval = 0;
    int node;

    cp_heur_env *heur_env = ea_env->heur;

    node = rand() % cp->n;

    if (cp->data->from == node || cp->data->to == node)
    {
        ;
    }
    else if (sol->selected[node])
    {
        cp_drop_sol_node(cp, heur_env, sol, node);
    }
    else
    {
        cp_add_sol_node(cp, heur_env, sol, node);
    }

    return rval;
}
