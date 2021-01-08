#include "../bac.h"
#include "op-solver.h"

static int
heuristic_clique_(solver_graph *graph, graph_vertex *u,
                  graph_clique_repo *repo);

int
cp_shrink_exact_bac_graph_s1(cp_prob *cp, cp_exact_bac_env *env,
                             solver_graph *graph, graph_vertex *qstart,
                             graph_clique_repo *repo)
{
    int rval = 0;
    graph_vertex *v, *u;
    graph_arc *e_vu = NULL, *h;
    graph_vertex *qhead, *qtail, *other;
    double c;
    int reorder;

    cp_exact_bac_param *bac_param = env->param;

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

        // par->count_queue++;

        c = v->y;

        u    = NULL;
        e_vu = NULL;
        for (e_vu = v->edge;
             e_vu && (fabs(e_vu->x - c) >= SOLVER_ZEROPLUS ||
                      fabs(e_vu->tail->y - e_vu->head->y) >= SOLVER_ZEROPLUS);
             e_vu = outnext(e_vu, v))
        {
            if (e_vu->x <= c - SOLVER_ZEROPLUS)
            {
                if (e_vu->x > e_vu->tail->y + SOLVER_ZEROPLUS)
                {
                    rval = heuristic_clique_(graph, e_vu->head, repo);
                }
                if (e_vu->x > e_vu->head->y + SOLVER_ZEROPLUS)
                {
                    rval = heuristic_clique_(graph, e_vu->tail, repo);
                }
            }
            if (bac_param->srk_s2 && graph->n3v > 2)
            {
                u = otherend(e_vu, v);
                if (e_vu->x > v->y + SOLVER_ZEROPLUS &&
                    e_vu->x > u->y + SOLVER_ZEROPLUS)
                {
                    if (v->i != graph->tail->i && u->i != graph->tail->i)
                        reorder = 0;
                    else
                        reorder = 1;

                    graph_identify_vertices(graph, v, u);
                    // bac_param->count_s2++;
                    ADD_TO_SRK_QUEUE(v);
                    for (h = v->edge; h; h = outnext(h, v))
                    {
                        other = otherend(h, v);
                        if (h->x > other->y + SOLVER_ZEROPLUS)
                        {
                            rval = clique_register_repo_srkvertices(
                            graph, v, other, h->x, repo);
                            check_rval(
                            rval, "clique_register_repo_srkvertices failed",
                            CLEANUP);
                        }
                        ADD_TO_SRK_QUEUE(other);
                    }

                    if (reorder)
                        graph_reorder_vertices(graph);

                    goto GET_OUT;
                }
            }
        }

        if (e_vu)
            u = otherend(e_vu, v);

        /* Rule S1 */
        if (e_vu && u->i != graph->tail->i)
        {
            graph_identify_vertices(graph, v, u);
            // par->count_s1++;

            ADD_TO_SRK_QUEUE(v);
            for (h = v->edge; h; h = outnext(h, v))
            {
                other = otherend(h, v);
                if (h->x > other->y + SOLVER_ZEROPLUS)
                {
                    rval = clique_register_repo_srkvertices(graph, v, other,
                                                            h->x, repo);
                    check_rval(rval, "clique_register_repo_srkvertices failed",
                               CLEANUP);
                }
                ADD_TO_SRK_QUEUE(other);
            }
        }
        else if (e_vu && u->i == graph->tail->i)
        { // To avoid reordering the nodes.
            ADD_TO_SRK_QUEUE(u);
        }
    GET_OUT:;
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

    // if (cvcount <= graph->nv / 2 && cvcount > 2 && cvcount < graph->nv - 2)
    if (cvcount <= graph->nv / 2)
        clique = clique_conv_vertices2clique(graph->orig, cverts, cvcount);
    else
        clique = clique_conv_vertices2coclique(graph->orig, cverts, cvcount);
    // clique->val = 2 * u->y + 2 * v->y - 2 * eweight;
    clique_register_repo(graph->orig, repo, clique);

CLEANUP:
    if (clique)
        clique_free(&clique);
    if (cverts)
        free(cverts);
    return rval;
}
