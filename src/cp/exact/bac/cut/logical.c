#include "cp/cp.h"
#include "ip/exact/bac/bac.h"
#include "cp/exact/bac/bac.h"

int
cp_sep_logical(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
               cp_cut **cuts)
{
    int rval = 0;
    int i;
    solver_graph *graph = bac_env->ip->lp->sol->graph;

    cp_cut_logic *logical;
    cp_cut *cut;

    *cutcount = 0;
    *cuts     = NULL;

    for (i = 0; i < graph->na; i++)
    {
        graph_arc *arc = graph->arcs[i];

        if (arc->x > arc->tail->y + SOLVER_IP_BAC_MIN_VIOL)
        {
            cut = cp_create_cut();
            check_null(cut, "", CLEANUP);
            cut->rhs   = 0.0;
            cut->sense = 'L';

            logical         = malloc(sizeof(cp_cut_logic));
            logical->arc[0] = arc->tail->i;
            logical->arc[1] = arc->head->i;
            logical->v      = arc->tail->i;

            cut->logical = logical;

            cut->next = *cuts;
            *cuts     = cut;
            (*cutcount)++;
        }
        if (arc->x > arc->head->y + SOLVER_IP_BAC_MIN_VIOL)
        {
            cut = cp_create_cut();
            check_null(cut, "", CLEANUP);
            cut->rhs   = 0.0;
            cut->sense = 'L';

            logical         = malloc(sizeof(cp_cut_logic));
            logical->arc[0] = arc->tail->i;
            logical->arc[1] = arc->head->i;
            logical->v      = arc->head->i;

            cut->logical = logical;

            cut->next = *cuts;
            *cuts     = cut;
            (*cutcount)++;
        }
    }

CLEANUP:
    return rval;
}
