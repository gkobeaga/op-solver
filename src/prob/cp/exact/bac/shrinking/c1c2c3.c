#include "../bac.h"
#include "op-solver.h"

int
cp_shrink_exact_bac_graph_c1c2c3(cp_prob *cp, cp_exact_bac_env *env,
                                 solver_graph *graph, graph_vertex *qstart,
                                 graph_clique_repo *repo)
{
    int rval = 0;
    graph_vertex *v, *u, *w, *t;
    graph_arc *e, *e_vu, *e_ut, *e_uw, *e_wt;
    graph_vertex *qhead, *qtail, *other;
    int try_c1c2, try_c3;
    double c;

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

        try_c1c2 = 0;
        try_c3   = 0;
        for (e_vu = v->edge; e_vu; e_vu = outnext(e_vu, v))
        {
            u       = otherend(e_vu, v);
            u->vt_x = e_vu->x;
            if ((!try_c1c2 || try_c3 < 2) && fabs(u->y - c) < SOLVER_ZEROPLUS)
            {
                if (!try_c1c2 && fabs(e_vu->x - c) < SOLVER_ZEROPLUS)
                    try_c1c2++;
                try_c3++;
            }
        }

        // Rule 1 and 2
        if (try_c1c2)
        {
            for (e_vu = v->edge; e_vu; e_vu = outnext(e_vu, v))
            {
                u = otherend(e_vu, v);
                if (fabs(u->y - c) < SOLVER_ZEROPLUS &&
                    fabs(e_vu->x - c) < SOLVER_ZEROPLUS)
                {
                    for (e_ut = u->edge; e_ut; e_ut = outnext(e_ut, u))
                    {
                        t = otherend(e_ut, u);
                        if (t->i != v->i && fabs(t->y - c) < SOLVER_ZEROPLUS)
                        { // x({u,v}:t)=c
                            if (fabs(e_ut->x + t->vt_x - c) < SOLVER_ZEROPLUS)
                            {
                                graph_identify_vertices(graph, v, u);
                                // if (t->vt_x < SOLVER_ZEROPLUS)
                                //    bac_param->count_c1++;
                                // else
                                //    bac_param->count_c2++;

                                ADD_TO_SRK_QUEUE(v);
                                for (e = v->edge; e; e = outnext(e, v))
                                {
                                    other = otherend(e, v);

                                    if (e->x > other->y + SOLVER_ZEROPLUS)
                                    {

                                        rval = clique_register_repo_srkvertices(
                                        graph, v, other, e->x, repo);
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

        /* Rule C3 */
        if (try_c3 > 1)
        {
            for (e_vu = v->edge; e_vu; e_vu = outnext(e_vu, v))
            {
                u = otherend(e_vu, v);
                if (fabs(u->y - c) < SOLVER_ZEROPLUS)
                {
                    for (e_ut = u->edge; e_ut; e_ut = outnext(e_ut, u))
                    {
                        t       = otherend(e_ut, u);
                        t->ut_x = e_ut->x;
                    }
                    for (e_uw = u->edge; e_uw; e_uw = outnext(e_uw, u))
                    {
                        w = otherend(e_uw, u);
                        if (w->i != v->i && fabs(w->y - c) < SOLVER_ZEROPLUS)
                        {
                            if (fabs(e_vu->x + e_uw->x + w->vt_x - 2 * c) <
                                SOLVER_ZEROPLUS)
                            {
                                for (e_wt = w->edge; e_wt;
                                     e_wt = outnext(e_wt, w))
                                {
                                    t = otherend(e_wt, w);
                                    if (t->i != v->i && t->i != u->i &&
                                        fabs(t->y - c) < SOLVER_ZEROPLUS)
                                    {
                                        if (fabs(e_wt->x + t->vt_x + t->ut_x -
                                                 c) < SOLVER_ZEROPLUS)
                                        {
                                            for (e_ut = u->edge; e_ut;
                                                 e_ut = outnext(e_ut, u))
                                            {
                                                t       = otherend(e_ut, u);
                                                t->ut_x = 0.0;
                                            }
                                            graph_identify_vertices(graph, v,
                                                                    u);
                                            graph_identify_vertices(graph, v,
                                                                    w);
                                            // bac_param->count_c3++;

                                            ADD_TO_SRK_QUEUE(v);
                                            for (e = v->edge; e;
                                                 e = outnext(e, v))
                                            {
                                                other = otherend(e, v);

                                                if (e->x >
                                                    other->y + SOLVER_ZEROPLUS)
                                                {
                                                    rval =
                                                    clique_register_repo_srkvertices(
                                                    graph, v, other, e->x,
                                                    repo);
                                                    check_rval(
                                                    rval,
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
                    for (e_ut = u->edge; e_ut; e_ut = outnext(e_ut, u))
                    {
                        t       = otherend(e_ut, u);
                        t->ut_x = 0.0;
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
