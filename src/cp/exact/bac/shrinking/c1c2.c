#include "../bac.h"
#include "op-solver.h"

static int
heuristic_clique_(solver_graph *graph, graph_vertex *u,
                  graph_clique_repo *repo);

int
cp_shrink_exact_bac_graph_c1c2(cp_prob *cp, cp_exact_bac_env *env,
                               solver_graph *graph, graph_vertex *qstart,
                               graph_clique_repo *repo)
{
    int rval = 0;
    graph_vertex *v, *u, *t;
    graph_arc *e, *e_vu, *e_ut, *h, *next;
    graph_vertex *qhead, *qtail, *other;
    int try_c1c2;
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

        // Cycle Rules
        try_c1c2 = 0;
        for (e_vu = v->edge; e_vu; e_vu = outnext(e_vu, v))
        {
            u       = otherend(e_vu, v);
            u->vt_x = e_vu->x;
            if (!try_c1c2 && fabs(u->y - c) < SOLVER_ZEROPLUS &&
                fabs(e_vu->x - c) < SOLVER_ZEROPLUS)
                try_c1c2++;
            if (e_vu->x <= c - SOLVER_ZEROPLUS)
            {
                if (e_vu->x > v->y + SOLVER_ZEROPLUS)
                {
                    rval = heuristic_clique_(graph, u, repo);
                }
                if (e_vu->x > u->y + SOLVER_ZEROPLUS)
                {
                    rval = heuristic_clique_(graph, v, repo);
                }
            }
        }

        if (try_c1c2)
        {
            for (e_vu = v->edge; e_vu; e_vu = outnext(e_vu, v))
            {
                u = otherend(e_vu, v);
                if (fabs(u->y - c) < SOLVER_ZEROPLUS &&
                    fabs(e_vu->x - c) < SOLVER_ZEROPLUS)
                { // Rule 1 and 2
                    for (e_ut = u->edge; e_ut; e_ut = next)
                    {
                        next = outnext(e_ut, u);
                        t    = otherend(e_ut, u);
                        if (t->i != v->i)
                        {
                            if (e_ut->x + t->vt_x >= c - SOLVER_ZEROPLUS)
                            {
                                if (e_ut->x > t->y + SOLVER_ZEROPLUS)
                                {
                                    rval = heuristic_clique_(graph, u, repo);
                                }
                                if (e_ut->x > u->y + SOLVER_ZEROPLUS)
                                {
                                    rval = heuristic_clique_(graph, t, repo);
                                }
                            }
                            if (t->i != v->i &&
                                fabs(t->y - c) < SOLVER_ZEROPLUS)
                            { // x({u,v}:t)=c

                                if (fabs(e_ut->x + t->vt_x - c) <
                                    SOLVER_ZEROPLUS)
                                {

                                    graph_identify_vertices(graph, v, u);
                                    // if (t->vt_x < SOLVER_ZEROPLUS)
                                    //    bac_param->count_c1++;
                                    // else
                                    //    bac_param->count_c2++;

                                    ADD_TO_SRK_QUEUE(v);
                                    for (h = v->edge; h; h = outnext(h, v))
                                    {
                                        other = otherend(h, v);

                                        if (h->x > other->y + SOLVER_ZEROPLUS)
                                        {
                                            rval =
                                            clique_register_repo_srkvertices(
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

static int
heuristic_clique_(solver_graph *graph, graph_vertex *u, graph_clique_repo *repo)
{
    int rval = 0;

    int cvcount;
    graph_vertex **cverts = NULL;
    graph_clique *clique  = NULL;

    rval = graph_expand_vertex(graph, u, &cvcount, &cverts);
    check_rval(rval, "graph_expand_vertex failed", CLEANUP);

    if (cvcount <= graph->nv / 2)
        clique = clique_conv_vertices2clique(graph->orig, cverts, cvcount);
    else
        clique = clique_conv_vertices2coclique(graph->orig, cverts, cvcount);
    clique_register_repo(graph->orig, repo, clique);

CLEANUP:
    if (clique)
        clique_free(&clique);
    if (cverts)
        free(cverts);
    return rval;
}
