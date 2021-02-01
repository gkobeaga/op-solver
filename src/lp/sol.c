#include "lp.h"

void
lp_create_sol_work(lp_sol *sol)
{
    sol->val      = 0.0;
    sol->x        = NULL;
    sol->rc       = NULL;
    sol->graph    = graph_create();
    sol->integral = 0;
    return;
}

lp_sol *
lp_create_sol(void)
{
    lp_sol *sol = malloc(sizeof(lp_sol));
    lp_create_sol_work(sol);
    return sol;
}

void
lp_free_sol(lp_sol **sol)
{
    if (*sol)
    {
        graph_free(&(*sol)->graph);
        free(*sol);
        *sol = NULL;
    }
    return;
}

int
lp_update_sol(lp_prob *lp, lp_sol *sol)
{
    int rval = 0;
    graph_arc *arc, *lparc, *next;
    graph_vertex *v     = NULL;
    graph_arc **xsorted = NULL;
    double *x           = NULL;
    graph_vertex **heap = NULL;
    int i;
    int narcs;

    sol->val = lp_get_objval(lp);

    check_assert(lp->graph->nv + lp->graph->na == lp_get_ncols(lp), "",
                 CLEANUP);

    x = malloc((lp->graph->nv + lp->graph->na) * sizeof(double));
    check_null(x, "out of memory\n", CLEANUP);

    rval = lp_get_x(lp, x);
    check_rval(rval, "failed\n", CLEANUP);

    graph_free(&sol->graph);
    sol->graph = graph_create();
    graph_add_vertices(sol->graph, lp->graph->nv);
    for (i = 0; i < lp->graph->nv; i++)
    {
        sol->graph->v[i]->fixed  = lp->graph->v[i]->fixed;
        sol->graph->v[i]->branch = lp->graph->v[i]->branch;
        sol->graph->v[i]->obj    = lp->graph->v[i]->obj;
    }
    sol->graph->orig  = lp->graph;
    sol->graph->fixed = NULL;
    for (v = lp->graph->fixed; v; v = v->fixed_next)
    {
        if (sol->graph->fixed)
            sol->graph->fixed->fixed_prev = sol->graph->v[v->i];
        sol->graph->v[v->i]->fixed_next = sol->graph->fixed;
        sol->graph->fixed               = sol->graph->v[v->i];
    }

    xsorted = malloc(lp->graph->na * sizeof(graph_arc *));
    check_null(xsorted, "out of memory\n", CLEANUP);

    sol->integral = 1;
    for (i = 0, narcs = 0; i < lp->graph->na; i++)
    {
        lparc = lp->graph->arcs[i];
        if (x[lp->graph->nv + lparc->ind] > SOLVER_ZEROPLUS)
        {
            xsorted[narcs]      = lp->graph->arcs[i];
            xsorted[narcs++]->x = x[sol->graph->nv + lparc->ind];
            if (x[lp->graph->nv + lparc->ind] <= SOLVER_ONEMINUS)
                sol->integral = 0;
        }
    }

    for (i = 0; i < narcs; i++)
    {
        lparc       = xsorted[i];
        arc         = graph_add_arc(sol->graph, lparc->tail->i, lparc->head->i);
        arc->ind    = lparc->ind;
        arc->cost   = lparc->cost;
        arc->branch = lparc->branch;
        arc->x      = x[sol->graph->nv + lparc->ind];
        arc->orig   = lparc;
        lparc->x    = 0.0;
        arc->tail->y += x[sol->graph->nv + lparc->ind] / 2;
        arc->head->y += x[sol->graph->nv + lparc->ind] / 2;
    }

    heap      = malloc(sol->graph->n3v * sizeof(graph_vertex *));
    int nheap = 0;
    for (v = sol->graph->tail; v; v = v->next)
    {
        v->active     = 1;
        heap[nheap++] = v;
    }
    while (nheap)
    {
        v         = heap[--nheap];
        v->active = 0;
        if (v->y > SOLVER_ZEROPLUS)
        {
            if (v->y <= SOLVER_ONEMINUS)
                sol->integral = 0;
        }
        else
        {
            if (v->deg)
            {
                for (arc = v->edge; arc; arc = next)
                {
                    graph_vertex *other = otherend(arc, v);
                    if (!other->active && other->deg > 1)
                    {
                        other->active = 1;
                        heap[nheap++] = other;
                    }
                    next = outnext(arc, v);
                    graph_del_arc(sol->graph, &arc);
                }
                sol->integral = 0;
            }
            v->y = 0;
            x[i] = 0;
        }
    }

#if 1
    if (sol->integral)
    {
        int ncomp;
        int *compscount      = NULL;
        graph_vertex **comps = NULL;
        rval = graph_get_comps(sol->graph, &ncomp, &compscount, &comps);
        check_rval(rval, "", CLEANUP);
        if (ncomp == 1)
            sol->graph->connected = 1;
        free(compscount);
        free(comps);
    }
#endif

    if (sol->x)
        free(sol->x);

    sol->x = x;

CLEANUP:
    if (rval && x)
        free(x);
    if (xsorted)
        free(xsorted);
    if (heap)
        free(heap);
    return rval;
}
