#include "../bac.h"
#include "op-solver.h"

int
cp_shrink_exact_bac_graph_c1(cp_prob *cp, cp_exact_bac_env *env,
                             solver_graph *graph, graph_vertex *qstart,
                             graph_clique_repo *repo)
{
    int rval = 0;
    graph_vertex *v, *u, *t;
    graph_arc *e, *e_vu, *e_ut, *h, *next;
    graph_vertex *qhead, *qtail, *other;
    int try_c1;
    double c;

    if (!graph->nv)
        return 0;

    if (qstart)
    {
        qhead = qstart;
        for (v = qstart; v->qnext; v = v->qnext) v->onqueue = 1;
        qtail          = v;
        qtail->onqueue = 1;
    }
    else
    {
        for (v = graph->tail; v->next; v = v->next)
        {
            v->qnext   = v->next;
            v->onqueue = 1;
        }
        qhead          = graph->tail;
        qtail          = v;
        qtail->onqueue = 1;
        qtail->qnext   = NULL;
    }

    while (qhead)
    {
        v     = qhead;
        qhead = qhead->qnext;
        if (!qhead)
            qtail = NULL;
        if (v->parent != v)
            continue;
        v->onqueue = 0;

        // bac_param->count_queue++;

        c = v->y;

        // Cycle Rules
        try_c1 = 0;
        for (e_vu = v->edge; e_vu; e_vu = outnext(e_vu, v))
        {
            u       = otherend(e_vu, v);
            u->vt_x = e_vu->x;
            if (!try_c1 && fabs(u->y - c) < SOLVER_ZEROPLUS &&
                fabs(e_vu->x - c) < SOLVER_ZEROPLUS)
                try_c1++;
        }

        if (try_c1)
        {
            for (e_vu = v->edge; e_vu; e_vu = outnext(e_vu, v))
            {
                u = otherend(e_vu, v);
                if (fabs(u->y - c) < SOLVER_ZEROPLUS &&
                    fabs(e_vu->x - c) < SOLVER_ZEROPLUS)
                { /* Rule C1 */
                    for (e_ut = u->edge; e_ut; e_ut = next)
                    {
                        next = outnext(e_ut, u);
                        t    = otherend(e_ut, u);
                        if (t->i != v->i && fabs(t->y - c) < SOLVER_ZEROPLUS &&
                            t->vt_x < SOLVER_ZEROPLUS)
                        {
                            if (fabs(e_ut->x - c) < SOLVER_ZEROPLUS)
                            {
                                graph_identify_vertices(graph, v, u);
                                // bac_param->count_c1++;

                                ADD_TO_SRK_QUEUE(v);
                                for (h = v->edge; h; h = outnext(h, v))
                                {
                                    other = otherend(h, v);
                                    if (h->x > other->y + SOLVER_ZEROPLUS)
                                    {
                                        rval = clique_register_repo_srkvertices(
                                        graph, v, other, h->x, repo);
                                        check_rval(rval,
                                                   "clique_register_repo_"
                                                   "srkvertices failed",
                                                   CLEANUP);
                                    }
                                    ADD_TO_SRK_QUEUE(other);
                                }
                                goto GET_OUT;
                            }
                        }
                    }
                }
            }
        }
    GET_OUT:
        for (e = v->edge; e; e = outnext(e, v))
        {
            other       = otherend(e, v);
            other->vt_x = 0.0;
        }
    }

CLEANUP:

    return rval;
}
