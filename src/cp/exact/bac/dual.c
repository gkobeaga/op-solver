#include "cp/cp.h"
#include "ip/exact/bac/bac.h"
#include "cp/exact/bac/bac.h"
#include <gmp.h>

#define BIG_PRICE_GEN 20000

typedef struct bigprlist
{
    graph_vertex *vert;
    graph_arc *edge;
    mpf_t rc;
    mpf_t inrc;
} bigprlist;

static int
pricing_duals(cp_prob *cp, cp_exact_bac_env *bac_env, mpf_t *node_pi,
              mpf_t *node_piest, mpf_t *node_picalc, double *d0_pi,
              mpf_t *cut_pi, mpf_t *node_rc, mpf_t *inside_rc, mpf_t *clique_pi,
              mpf_t rhs_sum),
price_verts(cp_prob *cp, mpf_t *node_rc, mpf_t penalty),
price_arcs(cp_prob *cp, solver_graph *graph, cp_cut_repo *cuts, int ecount,
           bigprlist **elist, mpf_t *node_pi, double d0_pi, mpf_t *clique_pi),
get_arcs_from_piest_graph(cp_prob *cp, cp_exact_bac_env *bac_env,
                          solver_graph *graph, int *ngen, bigprlist **genlist,
                          mpf_t *inside_rc, mpf_t *node_piest, double d0_pi,
                          int *last_key, int *finished),
get_arcs_from_piest_complete(cp_prob *cp, cp_exact_bac_env *bac_env,
                             solver_graph *graph, int *ngen,
                             bigprlist **genlist, mpf_t *inside_rc,
                             mpf_t *node_piest, double d0_pi, int *last_ikey,
                             int *last_jkey, int *finished);

int
cp_verify_branch_infeas(void *prob, void *env, int *yesno)
{
    int rval = 0;
    mpf_t exactbound, max, zero;

    cp_prob *cp               = (cp_prob *)prob;
    cp_exact_bac_env *bac_env = (cp_exact_bac_env *)env;

    mpf_init(max);
    mpf_init(zero);
    mpf_set_d(max, bac_env->ip->infeas_bound);

    *yesno = 0;
    mpf_init(exactbound);
    rval = cp_get_branch_dual_bound(cp, bac_env, exactbound);
    check_rval(rval, "failed", CLEANUP);

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        printf("Exact node bound: %f\n", mpf_get_d(exactbound));

    if (mpf_cmp(exactbound, max) > 0 || mpf_cmp(exactbound, zero) < 0)
    {
        printf("Problem is shown to be infeasible\n");
        *yesno                   = 1;
        bac_env->ip->upperboundN = SOLVER_MAXDOUBLE;
    }
    else
    {
        printf("Did not verify an infeasible LP\n");
        *yesno = 0;
    }

CLEANUP:

    mpf_clear(exactbound);
    mpf_clear(max);
    mpf_clear(zero);

    return rval;
}

int
cp_get_branch_dual_bound(void *prob, void *env, mpf_t bound)
{
    int prcount;
    bigprlist **prlist = NULL;
    mpf_t zero, penalty, rhs_sum, rc;
    mpf_t *node_pi     = NULL;
    mpf_t *node_piest  = NULL;
    mpf_t *node_picalc = NULL;
    mpf_t *clique_pi   = NULL;
    mpf_t *node_rc     = NULL;
    mpf_t *inside_rc   = NULL;
    mpf_t *cut_pi      = NULL;
    double d0_pi;
    int nbranch = 0;
    int i, j, finished;
    int rval = 0;
    int last_key, last_ikey, last_jkey;
    solver_graph *graph;
    graph_arc *tedge;

    cp_prob *cp               = (cp_prob *)prob;
    cp_exact_bac_env *bac_env = (cp_exact_bac_env *)env;

    lp_prob *lp       = bac_env->ip->lp;
    cp_cut_repo *cuts = bac_env->cuts;

    for (i = 0; i < bac_env->ip->history_depth; i++)
    {
        if (bac_env->ip->history[i]->edge)
            nbranch++;
    }

    prcount = 0;
    if (BIG_PRICE_GEN >= nbranch)
    {
        prlist = malloc(BIG_PRICE_GEN * sizeof(bigprlist *));
        for (i = 0; i < BIG_PRICE_GEN; i++)
        {
            prlist[i] = malloc(sizeof(bigprlist));
            mpf_init(prlist[i]->rc);
            mpf_init(prlist[i]->inrc);
        }
    }
    else
    {
        prlist = malloc(nbranch * sizeof(bigprlist *));
        for (i = 0; i < nbranch; i++)
        {
            prlist[i] = malloc(sizeof(bigprlist));
            mpf_init(prlist[i]->rc);
            mpf_init(prlist[i]->inrc);
        }
    }
    check_null(prlist, "out of memory", CLEANUP);

    graph = graph_create();
    graph_add_vertices(graph, cp->n);
    for (i = 0; i < graph->nv; i++) graph->v[i]->obj = cp->data->obj_node[i];

    node_pi = malloc(cp->n * sizeof(mpf_t));
    check_null(node_pi, "out of memory ", CLEANUP);
    node_piest = malloc(cp->n * sizeof(mpf_t));
    check_null(node_piest, "out of memory ", CLEANUP);
    node_rc = malloc(cp->n * sizeof(mpf_t));
    check_null(node_rc, "out of memory ", CLEANUP);
    node_picalc = malloc(cp->n * sizeof(mpf_t));
    check_null(node_picalc, "out of memory ", CLEANUP);
    inside_rc = malloc(cp->n * (cp->n - 1) / 2 * sizeof(mpf_t));
    check_null(inside_rc, "out of memory ", CLEANUP);
    for (i = 0; i < cp->n; i++)
    {
        mpf_init(node_pi[i]);
        mpf_init(node_piest[i]);
        mpf_init(node_picalc[i]);
        mpf_init(node_rc[i]);
        if (lp->status != SOLVER_LP_INFEASIBLE)
            mpf_set_d(node_rc[i], cp->data->obj_node[i]);
    }

    if (cuts->cliques->size)
    {
        clique_pi = malloc(cuts->cliques->size * sizeof(mpf_t));
        check_null(clique_pi, "out of memory ", CLEANUP);
        for (i = 0; i < cuts->cliques->size; i++) mpf_init(clique_pi[i]);
    }
    if (cuts->count)
    {
        cut_pi = malloc(cuts->count * sizeof(mpf_t));
        check_null(cut_pi, "out of memory", CLEANUP);
        for (i = 0; i < cuts->count; i++) mpf_init(cut_pi[i]);
        for (i = 0; i < cp->n * (cp->n - 1) / 2; i++) mpf_init(inside_rc[i]);
    }

    mpf_init(rhs_sum);
    rval = pricing_duals(cp, bac_env, node_pi, node_piest, node_picalc, &d0_pi,
                         cut_pi, node_rc, inside_rc, clique_pi, rhs_sum);
    check_rval(rval, "pricing_duals failed", CLEANUP);

    last_key  = 0;
    last_ikey = 0;
    last_jkey = 1;
    finished  = 0;
    mpf_init(zero);
    mpf_init(penalty);
    mpf_init(rc);

    price_verts(cp, node_rc, penalty);

    while (!finished)
    {
        if (cp->data->graph)
            rval = get_arcs_from_piest_graph(cp, bac_env, graph, &prcount,
                                             prlist, inside_rc, node_piest,
                                             d0_pi, &last_key, &finished);
        else
            rval = get_arcs_from_piest_complete(
            cp, bac_env, graph, &prcount, prlist, inside_rc, node_piest, d0_pi,
            &last_ikey, &last_jkey, &finished);
        check_rval(rval, "get_arcs_from_piest failed", CLEANUP);

        rval = price_arcs(cp, graph, cuts, prcount, prlist, node_picalc, d0_pi,
                          clique_pi);
        check_rval(rval, "price_arcs failed", CLEANUP);

        for (i = 0; i < prcount; i++)
        {
            if (mpf_cmp(prlist[i]->rc, zero) > 0)
            {
                mpf_add(penalty, penalty, prlist[i]->rc);
            }
            mpf_set(prlist[i]->rc, zero);
        }

#if 1
        graph_free(&graph);
        graph = graph_create();
        graph_add_vertices(graph, cp->n);
#else
        graph_delall_arcs(graph);
#endif
    }

    if (lp->status != SOLVER_LP_INFEASIBLE)
    {
        mpf_set(bound, zero);
        mpf_add(bound, bound, rhs_sum);
        mpf_add(bound, bound, penalty);

        if (cp->ip->upperboundG > mpf_get_d(bound) && mpf_get_d(bound) > 0)
            cp->ip->upperboundG = floor(mpf_get_d(bound));
        if (cp->ip->upperboundG < cp->ip->lowerboundG)
            cp->ip->upperboundG = cp->ip->lowerboundG;
    }

    if (nbranch)
    {
        graph_arc *edge;

        /* Adjust bound for branch edges                                  */

        prcount = 0;
#if 1
        graph_free(&graph);
        graph = graph_create();
        graph_add_vertices(graph, cp->n);
#else
        graph_delall_arcs(graph);
#endif
        for (i = 0; i < bac_env->ip->history_depth; i++)
        {
            if (bac_env->ip->history[i]->edge)
            {
                edge  = bac_env->ip->history[i]->edge;
                tedge = graph_find_arc(graph, edge->tail, edge->head);
                if (!tedge)
                {
                    tedge = graph_add_arc(graph, edge->tail->i, edge->head->i);
                    edge->aux = tedge->aux = prcount;

                    tedge->cost = (double)data_get_norm(
                    cp->data, tedge->tail->i, tedge->head->i);
                    mpf_set(prlist[prcount]->rc, zero);
                    int pos = inside_pos(cp->n, edge->tail->i, edge->head->i);
                    mpf_set(prlist[prcount]->inrc, inside_rc[pos]);
                    tedge->branch           = edge->branch;
                    prlist[prcount++]->edge = tedge;
                }
            }
        }

        rval = price_arcs(cp, graph, cuts, prcount, prlist, node_picalc, d0_pi,
                          clique_pi);
        check_rval(rval, "price_arcs failed", CLEANUP);

        for (j = 0; j < bac_env->ip->history_depth; j++)
        {
            while (!(bac_env->ip->history[j]->edge))
            {
                j++;
                if (j >= bac_env->ip->history_depth)
                    goto BRANCH_END;
            }

            edge = graph_find_arc(graph, bac_env->ip->history[j]->edge->tail,
                                  bac_env->ip->history[j]->edge->head);

            if (bac_env->ip->history[j]->rhs > 0)
            {
                mpf_set(rc, prlist[edge->aux]->rc);
                if (mpf_cmp(rc, zero) < 0)
                    mpf_add(penalty, penalty, rc);
            }
            else
            {
                mpf_set(rc, prlist[edge->aux]->rc);
                if (mpf_cmp(rc, zero) > 0)
                    mpf_sub(penalty, penalty, rc);
            }
        }
    }

BRANCH_END:

    mpf_set(bound, zero);
    mpf_add(bound, bound, rhs_sum);
    mpf_add(bound, bound, penalty);

    rval = 0;

CLEANUP:

    bac_env->ip->farkas_pricing = 0;

    if (BIG_PRICE_GEN >= nbranch)
    {
        for (i = 0; i < BIG_PRICE_GEN; i++)
        {
            mpf_clear(prlist[i]->rc);
            mpf_clear(prlist[i]->inrc);
            free(prlist[i]);
        }
    }
    else
    {
        for (i = 0; i < nbranch; i++)
        {
            mpf_clear(prlist[i]->rc);
            mpf_clear(prlist[i]->inrc);
            free(prlist[i]);
        }
    }

    graph_free(&graph);
    if (prlist)
        free(prlist);
    mpf_clear(zero);
    mpf_clear(rc);
    mpf_clear(penalty);
    mpf_clear(rhs_sum);
    for (i = 0; i < cuts->count; i++) mpf_clear(cut_pi[i]);
    if (cut_pi)
        free(cut_pi);
    for (i = 0; i < cuts->cliques->size; i++) mpf_clear(clique_pi[i]);
    if (clique_pi)
        free(clique_pi);
    for (i = 0; i < cp->n; i++)
    {
        mpf_clear(node_pi[i]);
        mpf_clear(node_piest[i]);
        mpf_clear(node_picalc[i]);
        mpf_clear(node_rc[i]);
    }
    if (node_pi)
        free(node_pi);
    if (node_rc)
        free(node_rc);
    if (node_piest)
        free(node_piest);
    if (node_picalc)
        free(node_picalc);
    if (cuts->count)
    {
        for (i = 0; i < cp->n * (cp->n - 1) / 2; i++) mpf_clear(inside_rc[i]);
    }
    if (inside_rc)
        free(inside_rc);

    return rval;
}

static int
get_pi(cp_prob *cp, cp_exact_bac_env *bac_env, mpf_t *node_pi, double *d0_pi,
       mpf_t *cut_pi)
{
    int i;
    int ncount        = cp->n;
    lp_prob *lp       = bac_env->ip->lp;
    cp_cut_repo *cuts = bac_env->cuts;
    int ncuts         = cuts->count;
    int nrows         = ncount + 1 + ncuts;
    double *pi        = (double *)NULL;
    int rval          = 0;

    pi = malloc(nrows * sizeof(double));
    check_null(pi, "out of memory ", CLEANUP);

    bac_env->ip->farkas_pricing = 0;
    rval                        = lp_get_pi(lp, pi);
    if (rval == 2)
    {
        bac_env->ip->farkas_pricing = 1;
        for (i = 0; i < nrows; i++) pi[i] *= -1;
        rval = 0;
    }
    check_rval(rval, "lp_get_pi failed", CLEANUP);

    for (i = 0; i < ncount; i++) mpf_set_d(node_pi[i], pi[i]);

    *d0_pi = pi[ncount];

    if (cut_pi)
    {
        for (i = 0; i < ncuts; i++) mpf_set_d(cut_pi[i], pi[ncount + 1 + i]);
    }

CLEANUP:
    free(pi);

    return rval;
}

static int
price_verts(cp_prob *cp, mpf_t *node_rc, mpf_t penalty)
{
    int rval = 0;
    int i;

    for (i = 0; i < cp->n; i++)
    {
        if (!i)
            mpf_add(penalty, penalty, node_rc[i]);
        else if (i && mpf_get_d(node_rc[i]) > 0)
            mpf_add(penalty, penalty, node_rc[i]);
    }

    return rval;
}

static int
pricing_duals(cp_prob *cp, cp_exact_bac_env *bac_env, mpf_t *node_pi,
              mpf_t *node_piest, mpf_t *node_picalc, double *d0_pi,
              mpf_t *cut_pi, mpf_t *node_rc, mpf_t *inside_rc, mpf_t *clique_pi,
              mpf_t rhs_sum)
{
    mpf_t x, zero, add, mul;
    int i, j, tmp;
    int rval = 0;
    int pos;

    lp_prob *lp       = bac_env->ip->lp;
    cp_cut_repo *cuts = bac_env->cuts;
    cp_cut *cut;

    mpf_init(zero);
    mpf_init(x);
    mpf_init(add);
    mpf_init(mul);

    rval = get_pi(cp, bac_env, node_pi, d0_pi, cut_pi);
    check_rval(rval, "get_pi failed", CLEANUP);

#if 0
  // RHS = 0
  for (i = 0; i < cp->n; i++)
  { mpf_set_d(add, 2*lp->result->graph->v[i]->x);
    mpf_mul (add, add, node_pi[i]);
    mpf_sub (rhs_sum, rhs_sum, add);
  }
#endif

    for (i = 0; i < cp->n; i++)
    {
        mpf_add(node_rc[i], node_rc[i], node_pi[i]);
        mpf_add(node_rc[i], node_rc[i], node_pi[i]);
    }

    // D0
    if (cp->data->cap)
    {
        mpf_set_d(x, *d0_pi);
        mpf_set_d(add, cp->data->cap);
        mpf_mul(add, add, x);
        mpf_add(rhs_sum, rhs_sum, add);
    }

    // CUTS
    for (i = 0; i < cuts->count; i++)
    {
        mpf_set(x, cut_pi[i]);
        mpf_set_d(add, cuts->cuts[i]->rhs);
        mpf_mul(add, add, x);
        mpf_add(rhs_sum, rhs_sum, add);
    }

    for (i = 0; i < cuts->cliques->size; i++) mpf_set_ui(clique_pi[i], 0);

    for (i = 0; i < cuts->count; i++)
    {
        cut = cuts->cuts[i];
        if (cut->tcount + cut->hcount)
        {
            for (j = 0; j < cut->hcount; j++)
            {
                mpf_sub(clique_pi[cut->handle_cid[j]],
                        clique_pi[cut->handle_cid[j]], cut_pi[i]);
            }

            for (j = 0; j < cut->tcount; j++)
            {
                mpf_sub(clique_pi[cut->teeth_cid[j]],
                        clique_pi[cut->teeth_cid[j]], cut_pi[i]);

                mpf_add(node_rc[cut->verts[2 * j]], node_rc[cut->verts[2 * j]],
                        cut_pi[i]);
                mpf_add(node_rc[cut->verts[2 * j]], node_rc[cut->verts[2 * j]],
                        cut_pi[i]);
                mpf_add(node_rc[cut->verts[2 * j + 1]],
                        node_rc[cut->verts[2 * j + 1]], cut_pi[i]);
                mpf_add(node_rc[cut->verts[2 * j + 1]],
                        node_rc[cut->verts[2 * j + 1]], cut_pi[i]);
            }
        }
        else if (cut->logical)
        {
            mpf_add(node_rc[cut->logical->v], node_rc[cut->logical->v],
                    cut_pi[i]);

            pos = inside_pos(cp->n, cut->logical->arc[0], cut->logical->arc[1]);
            mpf_sub(inside_rc[pos], inside_rc[pos], cut_pi[i]);
        }
        else if (cut->cover_edge)
        {
            for (j = 0; j < cut->cover_edge->na; j++)
            {
                pos = inside_pos(cp->n, cut->cover_edge->arcs[2 * j],
                                 cut->cover_edge->arcs[2 * j + 1]);
                mpf_sub(inside_rc[pos], inside_rc[pos], cut_pi[i]);
            }

            if (cut->cover_edge->strong)
            {
                lp->graph->marker++;
                for (j = 0; j < cut->cover_edge->na; j++)
                {
                    if (lp->graph->v[cut->cover_edge->arcs[2 * j]]->mark !=
                        lp->graph->marker)
                        mpf_add(node_rc[cut->cover_edge->arcs[2 * j]],
                                node_rc[cut->cover_edge->arcs[2 * j]],
                                cut_pi[i]);
                    if (lp->graph->v[cut->cover_edge->arcs[2 * j + 1]]->mark !=
                        lp->graph->marker)
                        mpf_add(node_rc[cut->cover_edge->arcs[2 * j + 1]],
                                node_rc[cut->cover_edge->arcs[2 * j + 1]],
                                cut_pi[i]);
                    lp->graph->v[cut->cover_edge->arcs[2 * j]]->mark =
                    lp->graph->marker;
                    lp->graph->v[cut->cover_edge->arcs[2 * j + 1]]->mark =
                    lp->graph->marker;
                }
            }
        }
        else if (cut->cover_vertex)
        {
            FOREACH_NODE_IN_CLIQUE (j, cut->cover_vertex->verts, tmp)
            {
                mpf_add(node_rc[j], node_rc[j], cut_pi[i]);
            }
        }
        else if (cut->path)
        {
            lp->graph->marker++;
            for (j = 0; j < cut->path->na; j++)
            {
                pos = inside_pos(cp->n, cut->path->arcs[2 * j],
                                 cut->path->arcs[2 * j + 1]);
                mpf_sub(inside_rc[pos], inside_rc[pos], cut_pi[i]);

                if (lp->graph->v[cut->path->arcs[2 * j]]->mark ==
                    lp->graph->marker)
                    mpf_add(node_rc[cut->path->arcs[2 * j]],
                            node_rc[cut->path->arcs[2 * j]], cut_pi[i]);
                if (lp->graph->v[cut->path->arcs[2 * j + 1]]->mark ==
                    lp->graph->marker)
                    mpf_add(node_rc[cut->path->arcs[2 * j + 1]],
                            node_rc[cut->path->arcs[2 * j + 1]], cut_pi[i]);
                lp->graph->v[cut->path->arcs[2 * j]]->mark = lp->graph->marker;
                lp->graph->v[cut->path->arcs[2 * j + 1]]->mark =
                lp->graph->marker;
            }
            for (j = 0; j < cut->path->fna; j++)
            {
                pos = inside_pos(cp->n, cut->path->farcs[2 * j],
                                 cut->path->farcs[2 * j + 1]);
                mpf_add(inside_rc[pos], inside_rc[pos], cut_pi[i]);
            }
        }
        else
        {
            printf("Unknown CUT %d\n", i);
            exit(1);
        }
    }

    for (i = 0; i < cp->n; i++)
    {
        mpf_set(node_piest[i], node_pi[i]);
        mpf_set(node_picalc[i], node_pi[i]);
    }

    for (i = 0; i < cuts->cliques->size; i++)
    {
        mpf_set(x, clique_pi[i]);
        if (cuts->cliques->cliques[i]->refcount)
        {
            if (mpf_cmp(x, zero) < 0)
            {
                FOREACH_NODE_IN_CLIQUE (j, cuts->cliques->cliques[i], tmp)
                {
                    mpf_sub(node_picalc[j], node_picalc[j], x);
                }
            }
            else if (mpf_cmp(x, zero) > 0)
            {
                FOREACH_NODE_IN_CLIQUE (j, cuts->cliques->cliques[i], tmp)
                {
                    mpf_sub(node_picalc[j], node_picalc[j], x);
                    mpf_sub(node_piest[j], node_piest[j], x);
                }
            }
        }
    }

CLEANUP:
    mpf_clear(add);
    mpf_clear(x);
    mpf_clear(zero);
    mpf_clear(mul);
    return rval;
}

static int
get_arcs_from_piest_graph(cp_prob *cp, cp_exact_bac_env *bac_env,
                          solver_graph *graph, int *ngen, bigprlist **genlist,
                          mpf_t *inside_rc, mpf_t *node_piest, double d0_pi,
                          int *last_key, int *finished)
{
    int i   = *last_key;
    int cnt = 0;
    double cost;
    mpf_t rc, zero, logrc;
    graph_arc *prarc, *lparc;

    *ngen     = 0;
    *finished = 0;
    mpf_init(rc);
    mpf_init(zero);
    mpf_init(logrc);

    if (!cp->data)
    {
        fprintf(stderr, "no source of edges in get_arcs_from_piest\n");
        return 1;
    }

    do
    {
        lparc = cp->data->graph->arcs[i];
        if (cp->data->cap)
        {
            cost =
            (double)data_get_norm(cp->data, lparc->tail->i, lparc->head->i);
            mpf_set_d(rc, -d0_pi * cost);
        }
        else
            mpf_set(rc, zero);
        mpf_sub(rc, rc, node_piest[lparc->tail->i]);
        mpf_sub(rc, rc, node_piest[lparc->head->i]);

        int pos;
        pos = inside_pos(cp->n, lparc->tail->i, lparc->head->i);

        mpf_add(rc, rc, inside_rc[pos]);
        if (mpf_cmp(rc, zero) >= 0)
        {
            prarc       = graph_add_arc(graph, lparc->tail->i, lparc->head->i);
            prarc->cost = cost;
            prarc->aux  = cnt;
            prarc->ind  = -1;
            mpf_set(genlist[cnt]->rc, zero);
            mpf_add(genlist[cnt]->rc, genlist[cnt]->rc, inside_rc[pos]);
            genlist[cnt++]->edge = prarc;
            if (cnt == BIG_PRICE_GEN)
                goto NOT_FINISHED;
        }
    } while ((i = (i + 1) % cp->n) != *last_key);

    *ngen     = cnt;
    *finished = 1;

    mpf_clear(logrc);
    mpf_clear(zero);
    mpf_clear(rc);
    return 0;

NOT_FINISHED:

    *finished = 0;
    *last_key = i;
    mpf_clear(rc);
    mpf_clear(logrc);
    return 0;
}

static int
get_arcs_from_piest_complete(cp_prob *cp, cp_exact_bac_env *bac_env,
                             solver_graph *graph, int *ngen,
                             bigprlist **genlist, mpf_t *inside_rc,
                             mpf_t *node_piest, double d0_pi, int *last_ikey,
                             int *last_jkey, int *finished)
{
    int i = *last_ikey;
    int j = *last_jkey;
    double cost;
    int cnt = 0;

    lp_prob *lp = bac_env->ip->lp;

    *ngen     = 0;
    *finished = 0;

    if (i >= cp->n)
    {
        i = 0;
        j = i + 1;
    }

    mpf_t rc, zero;
    graph_arc *arc, *prarc;

    mpf_init(rc);
    mpf_init(zero);

    for (; j < cp->n; j++)
    {
        if (cp->data->cap)
        {
            cost = (double)data_get_norm(cp->data, i, j);
            mpf_set_d(rc, -d0_pi * cost);
        }
        else
            mpf_set(rc, zero);
        mpf_sub(rc, rc, node_piest[i]);
        mpf_sub(rc, rc, node_piest[j]);

        int pos;
        pos = inside_pos(cp->n, i, j);
        mpf_add(rc, rc, inside_rc[pos]);

        if (mpf_cmp(rc, zero) >= 0)
        {
            prarc       = graph_add_arc(graph, i, j);
            prarc->cost = cost;
            prarc->aux  = cnt;
            arc         = graph_find_arc_hash(lp->graph->archash, i, j);
            if (arc)
                prarc->ind = arc->ind;
            else
                prarc->ind = -1;
            mpf_set(genlist[cnt]->rc, zero);
            mpf_set(genlist[cnt]->inrc, inside_rc[pos]);
            genlist[cnt++]->edge = prarc;
            if (cnt == BIG_PRICE_GEN)
                goto NOT_FINISHED;
        }
    }
    while ((i = (i + 1) % cp->n) != 0)
    {
        for (j = i + 1; j < cp->n; j++)
        {
            if (cp->data->cap)
            {
                cost = (double)data_get_norm(cp->data, i, j);
                mpf_set_d(rc, -d0_pi * cost);
            }
            else
                mpf_set(rc, zero);
            mpf_sub(rc, rc, node_piest[i]);
            mpf_sub(rc, rc, node_piest[j]);

            int pos;
            pos = inside_pos(cp->n, i, j);
            mpf_add(rc, rc, inside_rc[pos]);

            if (mpf_cmp(rc, zero) >= 0)
            {
                prarc       = graph_add_arc(graph, i, j);
                prarc->cost = cost;
                prarc->aux  = cnt;
                arc         = graph_find_arc_hash(lp->graph->archash, i, j);
                if (arc)
                    prarc->ind = arc->ind;
                else
                    prarc->ind = -1;
                mpf_set(genlist[cnt]->rc, zero);
                mpf_set(genlist[cnt]->inrc, inside_rc[pos]);
                genlist[cnt++]->edge = prarc;
                if (cnt == BIG_PRICE_GEN)
                    goto NOT_FINISHED;
            }
        }
    }

    *ngen     = cnt;
    *finished = 1;

    mpf_clear(zero);
    mpf_clear(rc);
    return 0;

NOT_FINISHED:

    *ngen      = cnt;
    *finished  = 0;
    *last_ikey = i;
    *last_jkey = j + 1;

    mpf_clear(zero);
    mpf_clear(rc);
    return 0;
}

static int
price_arcs(cp_prob *cp, solver_graph *graph, cp_cut_repo *cuts, int ecount,
           bigprlist **elist, mpf_t *node_pi, double d0_pi, mpf_t *clique_pi)
{
    int i, j, tmp;
    mpf_t x, zero, d0rhs;
    graph_vertex *v, *other;
    graph_arc *arc;
    graph_clique **cliques = cuts->cliques->cliques;
    int rval               = 0;
    bigprlist *bedge;

    mpf_init(zero);
    mpf_init(x);
    mpf_init(d0rhs);

    if (ecount == 0)
        goto CLEANUP;

    for (i = 0; i < ecount; i++)
    {
        bedge = elist[i];
        mpf_set(bedge->rc, zero);
        if (cp->data->cap)
        {
            mpf_set_d(d0rhs, (double)d0_pi * bedge->edge->cost);
            mpf_sub(bedge->rc, bedge->rc, d0rhs);
        }
        mpf_sub(bedge->rc, bedge->rc, node_pi[bedge->edge->tail->i]);
        mpf_sub(bedge->rc, bedge->rc, node_pi[bedge->edge->head->i]);
        mpf_add(bedge->rc, bedge->rc, bedge->inrc);
    }

    for (i = 0; i < cuts->cliques->size; i++)
    {
        mpf_set(x, clique_pi[i]);
        mpf_add(x, x, clique_pi[i]);
        if (mpf_cmp(x, zero) > 0)
        {
            graph->marker++;
            FOREACH_NODE_IN_CLIQUE (j, cliques[i], tmp)
            {
                v = graph->v[j];
                for (arc = v->edge; arc; arc = outnext(arc, v))
                {
                    other = otherend(arc, v);
                    if (other->mark == graph->marker && arc->aux != -1)
                        mpf_sub(elist[arc->aux]->rc, elist[arc->aux]->rc, x);
                }
                v->mark = graph->marker;
            }
        }
    }

CLEANUP:
    mpf_clear(x);
    mpf_clear(zero);
    mpf_clear(d0rhs);
    return rval;
}

int
cp_verify_branch_prune(void *prob, void *env, int *yesno)
{
    int rval = 0;
    mpf_t exactbound, lowerbound, one;

    cp_prob *cp               = (cp_prob *)prob;
    cp_exact_bac_env *bac_env = (cp_exact_bac_env *)env;

    mpf_init(exactbound);
    mpf_init(lowerbound);
    mpf_init(one);

    ip_prob *ip = cp->ip;

    *yesno = 0;
    rval   = cp_get_branch_dual_bound(cp, bac_env, exactbound);
    check_rval(rval, " failed", CLEANUP);

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        printf("Exact LP UPPER bound: %f\n", mpf_get_d(exactbound));

    mpf_set_d(lowerbound, ip->lowerboundG);

    mpf_set_ui(one, 1);
    mpf_add(lowerbound, lowerbound, one);

    if (mpf_cmp(exactbound, lowerbound) < 0)
    {
        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("Can prune lp.\n");

        *yesno = 1;
    }
    else
    {
        printf("Cannot prune lp.\n");
        printf("DUAL LP UPPER bound: %f\n", mpf_get_d(exactbound));
        printf("PRIMAL LP UPPER bound: %f\n", bac_env->ip->lp->sol->val);
        printf("Exact LP LOWER bound: %f\n", ip->lowerboundG);

        *yesno = 0;
    }

CLEANUP:

    mpf_clear(exactbound);
    mpf_clear(lowerbound);
    mpf_clear(one);

    return rval;
}
