#include "op-solver.h"

static void
add_primheap(solver_graph *graph, solver_dheap *dheap, int n, int *neighbor,
             int max, double sum_cost);

solver_graph *
graph_get_mst_max(solver_graph *graph)
{
    return graph_get_mst(graph, SOLVER_MST_MAX);
}

solver_graph *
graph_get_mst_min(solver_graph *graph)
{
    return graph_get_mst(graph, SOLVER_MST_MIN);
}

solver_graph *
graph_get_mst(solver_graph *graph, int max)
{
    int rval = 0;
    int i, n;
    int nvconn;
    int *neighbor   = (int *)NULL;
    double sum_cost = 0.0;
    solver_dheap *dheap;
    graph_arc *arc, *treearc;
    int comp = 0;

    solver_graph *tree = graph_create();
    graph_add_vertices(tree, graph->nv);

    neighbor = malloc(graph->nv * sizeof(int));
    check_null(neighbor, "Out of memory in solver_graph_spanningtree", CLEANUP);

    for (i = 1; i < graph->na; i++)
    {
        graph->arcs[i]->used = 0;
        if (max == SOLVER_MST_MAX)
            sum_cost += graph->arcs[i]->x;
    }

    for (i = 0, nvconn = 0; i < graph->nv; i++)
    {
        if (graph->v[i]->deg > 0)
            nvconn++;
        graph->v[i]->mark = graph->marker;
        tree->v[i]->y     = graph->v[i]->y;
        tree->v[i]->comp = graph->v[i]->comp = 0;
    }

    (graph->marker)++;

    dheap = dheap_create(graph->nv);
    for (i = 0; i < graph->nv; i++) neighbor[i] = -1;

    add_primheap(graph, dheap, 0, neighbor, max, sum_cost);
    graph->v[0]->mark = graph->marker;

    for (i = 1; i < nvconn; i++)
    {
        while (1)
        {
            n = dheap_deletemin(dheap);
            if (n == -1)
                goto CLEANUP;
            else if (neighbor[neighbor[n]] == -1)
                break;
            else
                add_primheap(graph, dheap, n, neighbor, max, sum_cost);
        }

        treearc = graph_add_arc(tree, n, neighbor[n]);
        arc     = graph_find_arc_hash(graph->archash, n, neighbor[n]);

        if (!n)
        {
            tree->v[neighbor[n]]->comp  = ++comp;
            graph->v[neighbor[n]]->comp = comp;
        }
        else
        {
            tree->v[neighbor[n]]->comp  = tree->v[n]->comp;
            graph->v[neighbor[n]]->comp = graph->v[n]->comp;
        }

        arc->used     = 1;
        treearc->ind  = arc->ind;
        treearc->cost = arc->cost;
        treearc->x    = arc->x;
        treearc->orig = arc;

        graph->v[neighbor[n]]->mark = graph->marker;

        add_primheap(graph, dheap, neighbor[n], neighbor, max, sum_cost);
        add_primheap(graph, dheap, n, neighbor, max, sum_cost);
    }

CLEANUP:
    dheap_free(&dheap);

    if (neighbor)
        free(neighbor);
    if (rval)
    {
        graph_free(&tree);
    }
    return tree;
}

static void
add_primheap(solver_graph *graph, solver_dheap *dheap, int n, int *neighbor,
             int max, double sum_cost)
{
    graph_vertex *v, *other;
    graph_arc *e;
    double minimum = SOLVER_MAXDOUBLE;

    v = graph->v[n];
    for (e = v->edge; e; e = outnext(e, v))
    {
        other = otherend(e, v);
        if (other->mark != graph->marker)
        {
            if (max == SOLVER_MST_MIN && e->x < minimum)
            {
                dheap->key[v->i] = e->x;
                minimum          = e->x;
                neighbor[v->i]   = other->i;
            }
            else if (max == SOLVER_MST_MAX && (sum_cost - e->x) < minimum)
            {
                dheap->key[v->i] = sum_cost - e->x;
                minimum          = sum_cost - e->x;
                neighbor[v->i]   = other->i;
            }
        }
    }

    if (minimum == SOLVER_MAXDOUBLE)
    {
        dheap->key[v->i] = SOLVER_MAXDOUBLE;
        neighbor[v->i]   = v->i;
    }
    else
        dheap_insert(dheap, n);
}
