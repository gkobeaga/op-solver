#include "op-solver.h"
#include "cp/cp.h"

int
cp_is_ipsol_integral_connected(ip_sol *sol)
{
    int rval = 0;
    int i, ncomp;
    graph_vertex **comps = NULL;
    int *compscount      = NULL;
    solver_graph *graph  = sol->graph;
    graph_vertex *v;
    graph_arc *arc;

    sol->integral    = 0;
    graph->connected = 0;
    for (i = 0; i < graph->nv; i++)
    {
        v = graph->v[i];
        if (v->y > 0.5)
        {
            if (v->y <= SOLVER_ONEMINUS)
                goto cleanup;
        }
        else
        {
            if (v->y > SOLVER_ZEROPLUS)
                goto cleanup;
        }
    }
    for (i = 0; i < graph->na; i++)
    {
        arc = graph->arcs[i];
        if (arc->x > 0.5)
        {
            if (arc->x <= SOLVER_ONEMINUS)
                goto cleanup;
        }
        else
        {
            if (arc->x > SOLVER_ZEROPLUS)
                goto cleanup;
        }
    }

    sol->integral = 1;

    graph_get_comps(graph, &ncomp, &compscount, &comps);

    if (ncomp > 1)
        goto cleanup;

    graph->connected = 1;
    rval             = 1;

cleanup:
    if (comps)
        free(comps);
    if (compscount)
        free(compscount);
    return rval;
}
