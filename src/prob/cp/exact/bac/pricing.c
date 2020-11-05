#include "op-solver.h"

static int
get_arcs_from_piest_graph(cp_prob *cp, cp_exact_bac_env *bac_env,
                          solver_graph *graph, int *narcs, graph_arc **arcs,
                          double d0_pi, double *node_piest, double *inside_rc,
                          int max_arcs, int *last_key, int *finished),
get_arcs_from_piest_complete(cp_prob *cp, cp_exact_bac_env *bac_env,
                             solver_graph *graph, int *narcs, graph_arc **arcs,
                             double d0_pi, double *node_piest,
                             double *inside_rc, int max_arcs, int start,
                             int *last_ikey, int *last_jkey, int *finished),
get_pi(cp_prob *cp, cp_exact_bac_env *bac_env, double *node_pi, double *d0_pi,
       double *cut_pi),
pricing_duals(cp_prob *cp, cp_exact_bac_env *bac_env, double *node_pi,
              double *node_piest, double *node_picalc, double *d0_pi,
              double *cut_pi, double *clique_pi, double *inside_rc);

static void
price_arcs(cp_prob *cp, solver_graph *graph, cp_cut_repo *cuts, int narcs,
           graph_arc **arcs, double *node_pi, double d0_pi, double *cut_pi,
           double *clique_pi, double *inside_rc);

static int
sort_bad(const void *xx, const void *yy)
/**************************************************************************/
{
    graph_arc *x = *(graph_arc **)xx, *y = *(graph_arc **)yy;

    if (x->rc < y->rc)
        return 1;
    if (x->rc > y->rc)
        return -1;

    return 0;
}

int
cp_recover_infeas(void *prob, void *env)
{
    int rval;
    int nadd, max_arcs;
    int finished;
    int last_key, last_ikey, last_jkey;
    graph_arc *arc;
    graph_arc **prlist  = NULL;
    graph_arc **estlist = NULL;
    int prcount;
    int estcount;
    double *node_pi     = NULL;
    double *node_piest  = NULL;
    double *node_picalc = NULL;
    double *clique_pi   = NULL;
    double *cut_pi      = NULL;
    double *inside_rc   = NULL;
    double d0_pi;
    int i, iend;
    int start                 = 0;
    cp_prob *cp               = (cp_prob *)prob;
    cp_exact_bac_env *bac_env = (cp_exact_bac_env *)env;
    solver_graph *graph;
    lp_prob *lp       = bac_env->ip->lp;
    cp_cut_repo *cuts = bac_env->cuts;

    double penalty;
    int nadded = 0;

    if (bac_env->ip->verbosity >= SOLVER_VERBOSITY_INFO)
        printf("infeas_recover ...\n");

    max_arcs =
    SOLVER_IP_BAC_PRICE_GEN + SOLVER_IP_BAC_PRICE_GEN_FACTOR * lp->graph->nv;

    stats_start(bac_env->stats->add_vars);

    prcount  = 0;
    estcount = 0;

    graph = graph_create();
    graph_add_vertices(graph, lp->graph->nv);

    inside_rc =
    malloc(lp->graph->nv * (lp->graph->nv - 1) / 2 * sizeof(double));
    memset(inside_rc, 0.0,
           lp->graph->nv * (lp->graph->nv - 1) / 2 * sizeof(double));

    estlist =
    malloc((SOLVER_IP_BAC_PRICE_POOL + max_arcs) * sizeof(graph_arc *));
    prlist =
    malloc((SOLVER_IP_BAC_PRICE_POOL + max_arcs) * sizeof(graph_arc *));
    node_pi     = malloc(lp->graph->nv * sizeof(double));
    node_piest  = malloc(lp->graph->nv * sizeof(double));
    node_picalc = malloc(lp->graph->nv * sizeof(double));
    if (!prlist || !node_pi || !node_piest)
    {
        fprintf(stderr, "out of memory \n");
        rval = 1;
        goto cleanup;
    }

    if (cuts->cliques->size)
    {
        clique_pi = malloc(cuts->cliques->size * sizeof(double));
        check_null(clique_pi, "out of memory", cleanup);
    }

    if (cuts->count)
    {
        cut_pi = malloc(cuts->count * sizeof(double));
        check_null(cut_pi, "out of memory", cleanup);
    }

    // Calc prices
    rval = pricing_duals(cp, bac_env, node_pi, node_piest, node_picalc, &d0_pi,
                         cut_pi, clique_pi, inside_rc);
    check_rval(rval, "pricing_duals failed.", cleanup);

    last_key  = 0;
    last_ikey = 0;
    last_jkey = 1;
    // Pricing LOOP
    finished = 0;
    nadded   = 0;
    penalty  = 0.0;

    while (!finished)
    {
        // Get a set of max_arcs arcs with rc_est >0
        if (0 && cp->data->graph)
            get_arcs_from_piest_graph(cp, bac_env, graph, &estcount, estlist,
                                      d0_pi, node_piest, inside_rc, max_arcs,
                                      &last_key, &finished);
        else
            get_arcs_from_piest_complete(
            cp, bac_env, graph, &estcount, estlist, d0_pi, node_piest,
            inside_rc, max_arcs, start, &last_ikey, &last_jkey, &finished);

        // Get prices
        price_arcs(cp, graph, cuts, estcount, estlist, node_picalc, d0_pi,
                   cut_pi, clique_pi, inside_rc);

        // Select the arcs whose price are positive
        for (i = 0; i < estcount; i++)
        {
            arc = estlist[i];
            if (arc->rc > 0.0)
                penalty += arc->rc;
            if (arc->rc > SOLVER_IP_BAC_PRICE_RCTHRESH)
                prlist[prcount++] = arc;
        }

        nadd = 0;
        while (
        (!finished && prcount >= SOLVER_IP_BAC_PRICE_POOL) ||
        (finished && penalty > SOLVER_IP_BAC_PRICE_MAXPENALTY && prcount > 0))
        {
            nadd = SOLVER_IP_PRICE_MAX_ADD;
            if (nadd >= prcount)
                nadd = prcount;
            else
                qsort(prlist, prcount, sizeof(graph_arc *), sort_bad);

            rval = cp_add_lp_arcs(cp, bac_env, prlist, nadd);
            check_rval(rval, " failed", cleanup);

            nadded += nadd;
            prcount -= nadd;

            // Dual optimize
            rval = lp_opt_dual(lp);
            if (lp->status != SOLVER_LP_INFEASIBLE)
            {
                if (bac_env->ip->verbosity >= SOLVER_VERBOSITY_INFO)
                    printf("LP is now feasible\n");
                rval = 1;
                goto DONE;
            }

            for (i = 0; i < nadd; i++)
            {
                graph_del_arc(graph, &(prlist[i]));
                prlist[i] = NULL;
            }

            for (i = nadd; i < prcount + nadd; i++)
            {
                prlist[i - nadd] = prlist[i];
                prlist[i]        = NULL;
            }

            // Calc duals
            rval = pricing_duals(cp, bac_env, node_pi, node_piest, node_picalc,
                                 &d0_pi, cut_pi, clique_pi, inside_rc);
            check_rval(rval, "pricing_duals failed.", cleanup);

            // Get prices
            price_arcs(cp, graph, cuts, prcount, prlist, node_picalc, d0_pi,
                       cut_pi, clique_pi, inside_rc);

            penalty = 0.0;
            for (i = 0, iend = prcount, prcount = 0; i < iend; i++)
            {
                arc = prlist[i];
                if (arc->rc > 0.0)
                    penalty += arc->rc;

                if (arc->rc > SOLVER_IP_BAC_PRICE_RCTHRESH)
                {
                    if (!graph_find_arc(lp->graph, arc->tail, arc->head))
                        prlist[prcount++] = arc;
                    else
                        graph_del_arc(graph, &arc);
                }
                else
                    graph_del_arc(graph, &arc);
            }
        }

        if (nadd > 0)
        {
            last_key  = 0;
            start     = last_ikey;
            last_jkey = start + 1;
            finished  = 0;
            estcount  = 0;
            prcount   = 0;
#if 1
            graph_free(&graph);
            graph = graph_create();
            graph_add_vertices(graph, lp->graph->nv);
#else
            graph_delall_arcs(graph);
#endif
            for (i = 0; i < prcount; i++)
                graph_add_arc(graph, prlist[i]->tail->i, prlist[i]->head->i);
        }
    }

DONE:

cleanup:

    if (bac_env->ip->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        if (lp->status != SOLVER_LP_INFEASIBLE)
            printf("Recovered a feasible LP\n");
        else
            printf("Could not recover a feasible LP\n");
    }

    rval = stats_stop(bac_env->stats->add_vars, nadded);

    if (cut_pi)
        free(cut_pi);
    if (clique_pi)
        free(clique_pi);
    if (node_piest)
        free(node_piest);
    if (node_pi)
        free(node_pi);
    if (node_picalc)
        free(node_picalc);
    if (estlist)
        free(estlist);
    if (prlist)
        free(prlist);
    if (inside_rc)
        free(inside_rc);

    graph_free(&graph);
    return rval;
}

int
cp_add_badvars(void *prob, void *env, int *nadded)
{
    int rval;
    int nadd, max_arcs;
    int finished;
    int last_key, last_ikey, last_jkey;
    graph_arc *arc;
    graph_arc **prlist  = NULL;
    graph_arc **estlist = NULL;
    int prcount, estcount;
    double *node_pi     = NULL;
    double *node_piest  = NULL;
    double *node_picalc = NULL;
    double *clique_pi   = NULL;
    double *cut_pi      = NULL;
    double *inside_rc   = NULL;
    double d0_pi;
    int i, iend;
    int start = 0;
    double penalty;

    cp_prob *cp               = (cp_prob *)prob;
    cp_exact_bac_env *bac_env = (cp_exact_bac_env *)env;
    solver_graph *graph;
    lp_prob *lp       = bac_env->ip->lp;
    cp_cut_repo *cuts = bac_env->cuts;

    *nadded = 0;

    max_arcs =
    SOLVER_IP_BAC_PRICE_GEN + SOLVER_IP_BAC_PRICE_GEN_FACTOR * lp->graph->nv;

    stats_start(bac_env->stats->add_vars);

    prcount  = 0;
    estcount = 0;

    graph = graph_create();
    graph_add_vertices(graph, lp->graph->nv);

    inside_rc =
    malloc(lp->graph->nv * (lp->graph->nv - 1) / 2 * sizeof(double));
    memset(inside_rc, 0.0,
           lp->graph->nv * (lp->graph->nv - 1) / 2 * sizeof(double));

    estlist =
    malloc((SOLVER_IP_BAC_PRICE_POOL + max_arcs) * sizeof(graph_arc *));
    prlist =
    malloc((SOLVER_IP_BAC_PRICE_POOL + max_arcs) * sizeof(graph_arc *));

    node_pi     = malloc(lp->graph->nv * sizeof(double));
    node_piest  = malloc(lp->graph->nv * sizeof(double));
    node_picalc = malloc(lp->graph->nv * sizeof(double));

    if (!prlist || !node_pi || !node_piest)
    {
        fprintf(stderr, "out of memory \n");
        rval = 1;
        goto cleanup;
    }

    if (cuts->cliques->size)
    {
        clique_pi = malloc(cuts->cliques->size * sizeof(double));
        check_null(clique_pi, "out of memory", cleanup);
    }

    if (cuts->count)
    {
        cut_pi = malloc(cuts->count * sizeof(double));
        check_null(cut_pi, "out of memory", cleanup);
    }

    // Calc prices
    rval = pricing_duals(cp, bac_env, node_pi, node_piest, node_picalc, &d0_pi,
                         cut_pi, clique_pi, inside_rc);
    check_rval(rval, "pricing_duals failed.", cleanup);

    last_key  = 0;
    last_ikey = 0;
    last_jkey = 1;
    // Pricing LOOP
    finished = 0;
    *nadded  = 0;
    penalty  = 0.0;

    while (!finished)
    { // Get a set of max_arcs arcs with rc_est >0
        if (cp->data->graph)
            get_arcs_from_piest_graph(cp, bac_env, graph, &estcount, estlist,
                                      d0_pi, node_piest, inside_rc, max_arcs,
                                      &last_key, &finished);
        else
            get_arcs_from_piest_complete(
            cp, bac_env, graph, &estcount, estlist, d0_pi, node_piest,
            inside_rc, max_arcs, start, &last_ikey, &last_jkey, &finished);

        // Get prices
        price_arcs(cp, graph, cuts, estcount, estlist, node_picalc, d0_pi,
                   cut_pi, clique_pi, inside_rc);

        // Select the arcs whose price are positive
        for (i = 0; i < estcount; i++)
        {
            arc = estlist[i];
            if (arc->rc > 0.0)
                penalty += arc->rc;
            if (arc->rc > SOLVER_IP_BAC_PRICE_RCTHRESH)
                prlist[prcount++] = arc;
        }

        nadd = 0;
        while (
        (!finished && prcount >= SOLVER_IP_BAC_PRICE_POOL) ||
        (finished && penalty > SOLVER_IP_BAC_PHASE1_MAXPENALTY && prcount > 0))
        {
            nadd = SOLVER_IP_PRICE_MAX_ADD;
            if (nadd >= prcount)
                nadd = prcount;
            else
                qsort(prlist, prcount, sizeof(graph_arc *), sort_bad);

            rval = cp_add_lp_arcs(cp, bac_env, prlist, nadd);
            check_rval(rval, " failed", cleanup);

            *nadded += nadd;
            prcount -= nadd;

            // Dual optimize
            rval = lp_opt_dual(lp);
            if (lp->status == SOLVER_LP_INFEASIBLE)
                printf("ADDBADS LP INFEASIBLE\n");
            if (lp->status == SOLVER_LP_INFEASIBLE)
            {
                fprintf(stderr, "Adding variables made LP infeasible!\n");
                rval = 1;
                goto cleanup;
            }
            else if (rval)
            {
                fprintf(stderr, "failed\n");
                goto cleanup;
            }

            for (i = 0; i < nadd; i++)
            {
                graph_del_arc(graph, &(prlist[i]));
                prlist[i] = NULL;
            }

            for (i = nadd; i < prcount + nadd; i++)
            {
                prlist[i - nadd] = prlist[i];
                prlist[i]        = NULL;
            }

            // Calc duals
            rval = pricing_duals(cp, bac_env, node_pi, node_piest, node_picalc,
                                 &d0_pi, cut_pi, clique_pi, inside_rc);
            check_rval(rval, "pricing_duals failed.", cleanup);

            // Get prices
            price_arcs(cp, graph, cuts, prcount, prlist, node_picalc, d0_pi,
                       cut_pi, clique_pi, inside_rc);

            penalty = 0.0;
            for (i = 0, iend = prcount, prcount = 0; i < iend; i++)
            {
                arc = prlist[i];
                if (arc->rc > 0.0)
                    penalty += arc->rc;

                if (arc->rc > SOLVER_IP_BAC_PRICE_RCTHRESH)
                {
                    if (!graph_find_arc(lp->graph, arc->tail, arc->head))
                        prlist[prcount++] = arc;
                    else
                        graph_del_arc(graph, &arc);
                }
                else
                    graph_del_arc(graph, &arc);
            }
        }

        if (nadd > 0)
        {
            last_key  = 0;
            start     = last_ikey;
            last_jkey = start + 1;
            finished  = 0;
            estcount  = 0;
            prcount   = 0;
#if 1
            graph_free(&graph);
            graph = graph_create();
            graph_add_vertices(graph, lp->graph->nv);
#else
            graph_delall_arcs(graph);
#endif
        }
    }

    if (*nadded)
        lp_update_sol(lp, lp->sol);

    double val;
    val = lp_get_objval(lp);
    check_assert(val >= 0, "lp_value failed", cleanup);

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        printf("(SPARSE) %d edges added, penalty %f, val %f\n", *nadded,
               penalty, val);

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        printf("New (node) upper bound: %f\n", val + penalty);
    bac_env->ip->upperboundN = val + penalty;

    if (bac_env->ip->depth == 0)
    {
        if (cp->ip->upperboundG > val + penalty)
            cp->ip->upperboundG = floor(val + penalty);
    }

cleanup:

    rval = stats_stop(bac_env->stats->add_vars, *nadded);

    if (cut_pi)
        free(cut_pi);
    if (clique_pi)
        free(clique_pi);
    if (node_piest)
        free(node_piest);
    if (node_pi)
        free(node_pi);
    if (node_picalc)
        free(node_picalc);
    if (estlist)
        free(estlist);
    if (prlist)
        free(prlist);
    if (inside_rc)
        free(inside_rc);

    graph_free(&graph);
    return rval;
}

static int
get_pi(cp_prob *cp, cp_exact_bac_env *bac_env, double *node_pi, double *d0_pi,
       double *cut_pi)
{
    int rval;
    double *pi        = NULL;
    lp_prob *lp       = bac_env->ip->lp;
    cp_cut_repo *cuts = bac_env->cuts;
    int nrows         = cp->n + 1 + cuts->count;
    int i;

    pi = malloc(nrows * sizeof(double));
    check_null(pi, "Out of memory in get_pi", cleanup);

    rval = lp_get_pi(lp, pi);
    if (rval == 2)
    {
        for (i = 0; i < nrows; i++) pi[i] *= -1;
        rval = 0;
    }
    check_rval(rval, "lp_get_pi failed", cleanup);

    if (node_pi)
    {
        for (i = 0; i < cp->n; i++) node_pi[i] = pi[i];
    }

    if (cp->data->cap)
        *d0_pi = pi[cp->n];

    if (cut_pi)
    {
        if (cp->data->cap)
            for (i = 0; i < cuts->count; i++) cut_pi[i] = pi[cp->n + 1 + i];
        else
            for (i = 0; i < cuts->count; i++) cut_pi[i] = pi[cp->n + i];
    }

    rval = 0;
cleanup:
    free(pi);
    return rval;
}

static int
pricing_duals(cp_prob *cp, cp_exact_bac_env *bac_env, double *node_pi,
              double *node_piest, double *node_picalc, double *d0_pi,
              double *cut_pi, double *clique_pi, double *inside_rc)
{
    double x;
    int i, tmp, j;
    int rval = 0;

    lp_prob *lp       = bac_env->ip->lp;
    cp_cut_repo *cuts = bac_env->cuts;

    cp_cut *cut;

    rval = get_pi(cp, bac_env, node_pi, d0_pi, cut_pi);
    check_rval(rval, "get_pi failed", cleanup);

    memset(clique_pi, 0.0, cuts->cliques->size * sizeof(double));

    memset(inside_rc, 0.0,
           lp->graph->nv * (lp->graph->nv - 1) / 2 * sizeof(double));

    for (i = 0; i < bac_env->cuts->count; i++)
    {
        cut = cuts->cuts[i];

        x = cut_pi[i];
        if (cut->tcount + cut->hcount)
        {
            for (j = 0; j < cut->hcount; j++)
                clique_pi[cut->handle_cid[j]] -= x;
            for (j = 0; j < cut->tcount; j++) clique_pi[cut->teeth_cid[j]] -= x;
        }
        else if ((cut->logical))
        {
            inside_rc[inside_pos(lp->graph->nv, cut->logical->arc[0],
                                 cut->logical->arc[1])] -= x;
        }
        else if ((cut->cover_edge))
        {
            for (j = 0; j < cut->cover_edge->na; j++)
                inside_rc[inside_pos(lp->graph->nv,
                                     cut->cover_edge->arcs[2 * j],
                                     cut->cover_edge->arcs[2 * j + 1])] -= x;
        }
        else if ((cut->path))
        {
            for (j = 0; j < cut->path->na; j++)
                inside_rc[inside_pos(lp->graph->nv, cut->path->arcs[2 * j],
                                     cut->path->arcs[2 * j + 1])] -= x;
            for (j = 0; j < cut->path->fna; j++)
                inside_rc[inside_pos(lp->graph->nv, cut->path->farcs[2 * j],
                                     cut->path->farcs[2 * j + 1])] += x;
        }
    }

    for (i = 0; i < lp->graph->nv; i++)
    {
        node_piest[i]  = node_pi[i];
        node_picalc[i] = node_pi[i];
    }

    for (i = 0; i < cuts->cliques->size; i++)
    {
        x = clique_pi[i];
        if (x < 0.0)
        {
            FOREACH_NODE_IN_CLIQUE (j, cuts->cliques->cliques[i], tmp)
            {
                node_picalc[j] -= x;
            }
        }
        else if (x > 0.0)
        {
            FOREACH_NODE_IN_CLIQUE (j, cuts->cliques->cliques[i], tmp)
            {
                node_picalc[j] -= x;
                node_piest[j] -= x;
            }
        }
    }

cleanup:
    return rval;
}

static int
get_arcs_from_piest_graph(cp_prob *cp, cp_exact_bac_env *bac_env,
                          solver_graph *graph, int *narcs, graph_arc **arcs,
                          double d0_pi, double *node_piest, double *inside_rc,
                          int max_arcs, int *last_key, int *finished)
{
    graph_arc *lparc, *prarc;
    int i = *last_key;
    double cost, rc_est;
    int cnt     = 0;
    lp_prob *lp = bac_env->ip->lp;

    *narcs    = 0;
    *finished = 0;

    if (!cp->data->graph->arcs)
        fprintf(stderr, "no source of edges to generate\n");

    // Full pricing
    do
    {
        lparc = cp->data->graph->arcs[i];
        if (cp->data->cap)
        {
            cost =
            (double)data_get_norm(cp->data, lparc->tail->i, lparc->head->i);
            rc_est = -d0_pi * cost;
        }
        else
            rc_est = 0;
        rc_est += -node_piest[lparc->tail->i] - node_piest[lparc->head->i];
        rc_est +=
        inside_rc[inside_pos(graph->nv, lparc->tail->i, lparc->head->i)];

        if (rc_est > 0.0)
        {
            if (!graph_find_arc(lp->graph, lparc->tail, lparc->head))
            {
                prarc = graph_add_arc(graph, lparc->tail->i, lparc->head->i);
                prarc->cost = cost;
                arcs[cnt++] = prarc;
                if (cnt == max_arcs)
                    goto NOT_FINISHED;
            }
        }
    } while ((i = (i + 1) % cp->n) != *last_key);

    *finished = 1;
    *narcs    = cnt;
    return 0;

NOT_FINISHED:

    *narcs    = cnt;
    *finished = 0;
    *last_key = i;
    return 0;
}

static int
get_arcs_from_piest_complete(cp_prob *cp, cp_exact_bac_env *bac_env,
                             solver_graph *graph, int *narcs, graph_arc **arcs,
                             double d0_pi, double *node_piest,
                             double *inside_rc, int max_arcs, int start,
                             int *last_ikey, int *last_jkey, int *finished)
{
    graph_arc *prarc;
    int i = *last_ikey;
    int j = *last_jkey;
    double cost, rc_est;
    int cnt = 0;

    lp_prob *lp = bac_env->ip->lp;

    *narcs    = 0;
    *finished = 0;

    if (i >= cp->n)
    {
        i = 0;
        j = i + 1;
    }

    for (; j < cp->n; j++)
    {

        if (cp->data->cap)
        {
            cost   = (double)data_get_norm(cp->data, i, j);
            rc_est = -d0_pi * cost;
        }
        else
            rc_est = 0;
        rc_est += -node_piest[i] - node_piest[j];
        rc_est += inside_rc[inside_pos(cp->n, i, j)];

        if (rc_est > 0.0)
        {
            if (!graph_find_arc_hash(lp->graph->archash, i, j))
            {
                prarc       = graph_add_arc(graph, i, j);
                prarc->cost = cost;
                arcs[cnt++] = prarc;
                if (cnt == max_arcs)
                    goto NOT_FINISHED;
            }
        }
    }
    while ((i = (i + 1) % cp->n) != start)
    {
        for (j = i + 1; j < cp->n; j++)
        {

            if (cp->data->cap)
            {
                cost   = (double)data_get_norm(cp->data, i, j);
                rc_est = -d0_pi * cost;
            }
            else
                rc_est = 0;
            rc_est += -node_piest[i] - node_piest[j];
            rc_est += inside_rc[inside_pos(cp->n, i, j)];

            if (rc_est > 0.0)
            {
                if (!graph_find_arc_hash(lp->graph->archash, i, j))
                {
                    prarc       = graph_add_arc(graph, i, j);
                    prarc->cost = cost;
                    arcs[cnt++] = prarc;
                    if (cnt == max_arcs)
                        goto NOT_FINISHED;
                }
            }
        }
    }

    *finished = 1;
    *narcs    = cnt;
    return 0;

NOT_FINISHED:

    *narcs     = cnt;
    *finished  = 0;
    *last_ikey = i;
    *last_jkey = j + 1;
    return 0;
}

static void
price_arcs(cp_prob *cp, solver_graph *graph, cp_cut_repo *cuts, int narcs,
           graph_arc **arcs, double *node_pi, double d0_pi, double *cut_pi,
           double *clique_pi, double *inside_rc)
{
    int i, j, tmp;
    double x;
    graph_vertex *v, *other;
    graph_arc *arc;

    if (narcs == 0)
        return;

    for (i = 0; i < narcs; i++)
    {
        arc = arcs[i];

        if (cp->data->cap)
        {
            arc->rc = -d0_pi * arc->cost;
        }
        else
            arc->rc = 0;
        arc->rc += -node_pi[arc->tail->i] - node_pi[arc->head->i];
        arc->rc += inside_rc[inside_pos(graph->nv, arc->tail->i, arc->head->i)];
    }

    for (i = 0; i < cuts->cliques->size; i++)
    {
        if (clique_pi[i] > 0)
        {
            x = clique_pi[i] * 2;
            graph->marker++;
            FOREACH_NODE_IN_CLIQUE (j, cuts->cliques->cliques[i], tmp)
            {
                v = graph->v[j];
                for (arc = v->edge; arc; arc = outnext(arc, v))
                {
                    other = otherend(arc, v);
                    if (other->mark == graph->marker)
                        arc->rc -= x;
                }
                v->mark = graph->marker;
            }
        }
    }
}
