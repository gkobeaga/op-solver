#include "lp/lp.h"
#include "op-solver.h"

static void
lp_create_prob_work(lp_prob *lp)
{
    lp->graph      = graph_create();
    lp->solver     = lp_create_solver();
    lp->infeasible = 0;
    lp->status     = SOLVER_LP_UNKNOWN;
    lp->sol        = lp_create_sol();
}

lp_prob *
lp_create_prob(void)
{
    lp_prob *lp = malloc(sizeof(lp_prob));
    lp_create_prob_work(lp);
    return lp;
}

void
lp_free_prob(lp_prob **lp)
{

    if (!(*lp))
        return;

    lp_free_solver(&(*lp)->solver);

    graph_free(&(*lp)->sol->graph);
    free((*lp)->sol->x);
    free((*lp)->sol);

    graph_free(&(*lp)->graph);

    free(*lp);
}

int
lp_set_graph(lp_prob *lp, solver_graph *graph)
{
    int rval = 0;
    int i;
    graph_arc *arc, *lparc;

    graph_add_vertices(lp->graph, graph->nv);

    for (i = 0; i < graph->nv; i++)
    {
        lp->graph->v[i]->obj = graph->v[i]->obj;
    }

    for (i = 0; i < graph->na; i++)
    {
        arc         = graph->arcs[i];
        lparc       = graph_add_arc(lp->graph, arc->tail->i, arc->head->i);
        lparc->cost = arc->cost;
        lparc->obj  = arc->obj;
        lp->graph->arcs[i] = lparc;
    }

    return rval;
}
