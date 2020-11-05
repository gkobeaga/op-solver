#include "op-solver.h"

static void
connect_search(solver_graph *graph, int n, int comp, int *dstack);

int
graph_get_comps(solver_graph *graph, int *ncomp, int **compscount,
                graph_vertex ***comps)
{
    int i, k, rval;
    int *nverts = NULL;
    int *dstack = NULL;
    graph_vertex *v;

    *comps      = NULL;
    *compscount = NULL;
    *ncomp      = 0;

    for (i = 0; i < graph->nv; i++) graph->v[i]->comp = 0;

    *comps = malloc(graph->nv * sizeof(graph_vertex *));
    check_null(*comps, "out of memory ", CLEANUP);

    dstack = malloc(graph->nv * sizeof(int));
    check_null(dstack, "out of memory ", CLEANUP);

    for (i = 0; i < graph->nv; i++) graph->v[i]->comp = 0;

    for (i = 0; i < graph->nv; i++)
    {
        if (!graph->v[i]->comp && graph->v[i]->y > SOLVER_ZEROPLUS)
        {
            if (graph->v[i]->edge)
            {
                (*ncomp)++;
                connect_search(graph, i, *ncomp, dstack);
            }
        }
    }

    *compscount = malloc((*ncomp) * sizeof(int));
    check_null(*compscount, "out of memory", CLEANUP);
    nverts = malloc((*ncomp) * sizeof(int));
    check_null(nverts, "out of memory ", CLEANUP);

    for (i = 0; i < *ncomp; i++) nverts[i] = 0;

    graph->marker++;
    for (i = 0; i < graph->nv; i++)
    {
        if (graph->v[i]->comp)
        {
            graph->v[i]->mark = graph->marker;
            nverts[graph->v[i]->comp - 1]++;
            for (v = graph->v[i]->members; v; v = v->members)
            {
                nverts[graph->v[i]->comp - 1]++;
            }
        }
    }

    for (i = 0, k = 0; i < *ncomp; i++)
    {
        (*compscount)[i] = nverts[i];
        nverts[i]        = k;
        k += (*compscount)[i];
    }

    for (i = 0; i < graph->nv; i++)
    {
        if (graph->v[i]->mark == graph->marker)
        {
            (*comps)[(nverts[graph->v[i]->comp - 1])++] = graph->v[i];
            for (v = graph->v[i]->members; v; v = v->members)
            {
                v->comp = graph->v[i]->comp;

                (*comps)[(nverts[graph->v[i]->comp - 1])++] = v;
            }
        }
    }

    rval = 0;

CLEANUP:

    if (nverts)
        free(nverts);
    if (dstack)
        free(dstack);

    if (rval && *compscount)
        free(*compscount);
    if (rval && *comps)
        free(*comps);

    return rval;
}

static void
connect_search(solver_graph *graph, int n, int comp, int *dstack)
{
    int head = 0;
    graph_vertex *v, *other;
    graph_arc *arc;

    graph->v[n]->comp = comp;
    dstack[head++]    = n;

    while (head > 0)
    {
        v = graph->v[dstack[--head]];
        for (arc = v->edge; arc; arc = outnext(arc, v))
        {
            other = otherend(arc, v);
            if (!other->comp && other->y > SOLVER_ZEROPLUS)
            {
                other->comp    = comp;
                dstack[head++] = other->i;
            }
        }
    }
}
