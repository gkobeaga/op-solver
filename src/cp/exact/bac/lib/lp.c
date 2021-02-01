#include "cp/cp.h"
#include "ip/exact/bac/bac.h"
#include "cp/exact/bac/bac.h"

static int
process_cuts(cp_prob *cp, cp_exact_bac_env *bac_env, int *nadded);

int
cp_add_lp_arcs(cp_prob *cp, cp_exact_bac_env *bac_env, graph_arc **arcs,
               int narcs)
{
    int i;
    int rval, nadd, oldnarc;
    graph_arc *arc;

    lp_prob *lp = bac_env->ip->lp;

    oldnarc = lp->graph->na;
    for (i = 0; i < narcs; i++)
    {
        arc = graph_add_arc(lp->graph, arcs[i]->tail->i, arcs[i]->head->i);
        arc->cost = arcs[i]->cost;
        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("      : adding arc (%d, %d) with index %d to the LP\n",
                   arc->tail->i, arc->head->i,

                   lp->graph->nv + arc->ind);
        }
    }
    nadd = lp->graph->na - oldnarc;

    lp_data *data;

    if (nadd)
    {
        data = cp_build_lp_data(cp, bac_env, lp->graph->nv + oldnarc,
                                lp->graph->nv + lp->graph->na - 1);
        check_null(data, " failed", cleanup);

        rval = lp_add_cols(lp, data);
        if (rval)
            goto cleanup;
    }
    else
        fprintf(stderr, "WARNING: none arcs add to LP\n");

    rval = 0;

cleanup:
    if (nadd)
    {
        if (data)
        {
            lp_free_data(&data);
        }
    }
    return rval;
}

lp_data *
cp_build_lp_data(cp_prob *cp, cp_exact_bac_env *bac_env, int col_start,
                 int col_end)
{
    int rval = 0;
    int i, j, ind;
    graph_vertex *vert;
    graph_arc *arc;
    graph_arc **nzlist;
    lp_data *data;
    int nzcnt;
    int nadd, ntotadd;

    lp_prob *lp         = bac_env->ip->lp;
    solver_graph *graph = lp->graph;
    cp_cut_repo *cuts   = bac_env->cuts;

    nadd = ntotadd = col_end - col_start + 1;

    if (col_start == 0)
        data = lp_create_data(graph->na + graph->nv, graph->nv + 1, 0);
    else
        data = lp_create_data(nadd, 0, 0);
    check_null(data, "out of memory", cleanup);

    // 5.1.6 Nodes
    if (col_start == 0)
    {
        for (i = 0; i < graph->nv; i++)
        {
            vert         = graph->v[i];
            data->obj[i] = vert->obj;

            if (vert->fixed > 0 || vert->branch > 0)
                data->lb[i] = 1.0;
            else
                data->lb[i] = 0.0;

            if (vert->fixed < 0 || vert->branch < 0)
                data->ub[i] = 0.0;
            else
                data->ub[i] = 1.0;

            data->cnt[i] = 1;
        }
        ntotadd = nadd;
        nadd    = nadd - graph->nv;

        // 5.1.5 Degree equations
        for (i = 0; i < graph->nv; i++)
        {
            data->rhs[i]   = 0.0;
            data->sense[i] = 'E';
        }
    }
    else if (col_start >= graph->nv)
    {
        col_start -= graph->nv;
    }
    else
    {
        fprintf(stderr, "col_start %d must be 0 or greater than graph->nv %d\n",
                col_start, graph->nv);
    }

    // 5.1.6 Arcs
    for (i = col_start; i < col_start + nadd; i++)
    {
        arc = graph->arcs[i];
        ind = (col_end - ntotadd + 1) < graph->nv ? graph->nv + arc->ind
                                                  : arc->ind - col_start;
        data->obj[ind] = 0.0;

        if (arc->fixed > 0 || arc->branch > 0)
            data->lb[ind] = 1.0;
        else
            data->lb[ind] = 0.0;

        if (arc->fixed < 0 || arc->branch < 0)
            data->ub[ind] = 0.0;
        else
            data->ub[ind] = 1.0;

        data->cnt[ind] = 3;
    }

    if (cuts->count)
    {
        nzlist = malloc(graph->na * sizeof(graph_arc *));
        for (i = 0; i < cuts->count; i++)
        {
            nzcnt = 0;
            cp_get_cut_arcs(graph, cuts->cuts[i], &nzcnt, nzlist);

            for (j = 0; j < nzcnt; j++)
            {
                arc = nzlist[j];
                ind = (col_end - ntotadd + 1) < graph->nv
                      ? graph->nv + arc->ind
                      : arc->ind - col_start;
                if (arc->coef != 0 && ind >= 0)
                    data->cnt[ind]++;
                arc->coef = 0;
            }
        }
        free(nzlist);
    }

    nzcnt = 0;
    for (i = 0; i < ntotadd; i++)
    {
        data->beg[i] = nzcnt;
        nzcnt += data->cnt[i];
        data->cnt[i] = 0;
    }

    data->ind = malloc(nzcnt * sizeof(int));
    data->val = malloc(nzcnt * sizeof(double));
    if (!data->ind || !data->val)
    {
        rval = 1;
        goto cleanup;
    }

    if (!(col_end - ntotadd + 1))
    {
        for (i = 0; i < graph->nv; i++)
        {
            graph_vertex *vert            = graph->v[i];
            data->val[data->beg[vert->i]] = -2.0;
            data->ind[data->beg[vert->i]] = vert->i;
            data->cnt[i]++;
        }
    }

    data->nzcnt = 0;
    for (i = col_start; i < col_start + nadd; i++)
    {
        arc = graph->arcs[i];
        ind = (col_end - ntotadd + 1) < graph->nv ? graph->nv + arc->ind
                                                  : arc->ind - col_start;
        data->val[data->beg[ind] + data->cnt[ind]] = 1.0;
        data->ind[data->beg[ind] + data->cnt[ind]] = arc->tail->i;
        data->cnt[ind]++;

        data->val[data->beg[ind] + data->cnt[ind]] = 1.0;
        data->ind[data->beg[ind] + data->cnt[ind]] = arc->head->i;
        data->cnt[ind]++;

        data->val[data->beg[ind] + data->cnt[ind]] =
        (double)data_get_norm(cp->data, arc->tail->i, arc->head->i);
        data->ind[data->beg[ind] + data->cnt[ind]] = graph->nv;
        data->cnt[ind]++;
        data->nzcnt += 3;
    }

    if (cuts->count)
    {
        nzlist = malloc(graph->na * sizeof(graph_arc *));
        for (i = 0; i < cuts->count; i++)
        {
            nzcnt = 0;
            cp_get_cut_arcs(graph, cuts->cuts[i], &nzcnt, nzlist);

            for (j = 0; j < nzcnt; j++)
            {
                arc = nzlist[j];
                ind = (col_end - ntotadd + 1) < graph->nv
                      ? graph->nv + arc->ind
                      : arc->ind - col_start;
                if (arc->coef != 0 && ind >= 0)
                {
                    data->val[data->beg[ind] + data->cnt[ind]] = arc->coef;
                    data->ind[data->beg[ind] + data->cnt[ind]] =
                    graph->nv + 1 + i;
                    data->cnt[ind]++;
                    data->nzcnt++;
                }
                arc->coef = 0;
            }
        }
        free(nzlist);
    }

    data->objsense = -1;
    if ((col_end - ntotadd + 1) == 0)
    {
        data->ncols = graph->na + graph->nv;
        data->nrows = graph->nv + 1 + cuts->count;
    }
    else
    {
        data->ncols = nadd;
        data->nrows = 0;
    }

    if (data->ncols > SOLVER_IP_MAX_COLS)
    {
        fprintf(stderr, " too much columns ");
        rval = 1;
        goto cleanup;
    }

    return data;

cleanup:

    if (rval)
        lp_free_data(&data);

    return NULL;
}

int
cp_add_lp_cuts(cp_prob *cp, cp_exact_bac_env *bac_env, cp_cut **cuts,
               int *nadded)
{
    int rval = 0;
    double objval;

    lp_prob *lp = bac_env->ip->lp;

    *nadded = 0;

    if (bac_env->cut_queue)
        bac_env->cut_queue->next = *cuts;
    else
        bac_env->cut_queue = *cuts;

    *cuts = NULL;

    rval = stats_start(bac_env->stats->add_cuts);
    check_rval(rval, "failed", cleanup);
    rval = process_cuts(cp, bac_env, nadded);
    check_rval(rval, "failed", cleanup);
    rval = stats_stop(bac_env->stats->add_cuts, *nadded);
    check_rval(rval, "failed", cleanup);

    if (bac_env->ip->lp->status == SOLVER_LP_SUCCESS)
    {
        objval = lp_get_objval(lp);
        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("  Add %2d cuts (Total %d), LP: %f (%.2ld seconds)\n",
                   *nadded, bac_env->cuts->count, objval,
                   stats_get_current_time(bac_env->stats->add_cuts));
        }
    }
    else if (bac_env->ip->lp->status == SOLVER_LP_INFEASIBLE)
    {
        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("  Problem infeasible after the cuts\n");
    }

cleanup:

    return rval;
}

static int
process_cuts(cp_prob *cp, cp_exact_bac_env *bac_env, int *nadded)
{
    int tnadded  = 0;
    int addpiece = 0;
    int rval;
    cp_cut *cut;
    lp_data *data = NULL;
    lp_prob *lp   = bac_env->ip->lp;
    double oldval = lp->sol->val;

    *nadded = 0;

    tnadded  = 0;
    addpiece = 0;
    while (bac_env->cut_queue)
    {
        cut                = bac_env->cut_queue;
        bac_env->cut_queue = cut->next;

        if (cp_eval_cut(lp->sol->graph, cut) < -SOLVER_IP_BAC_MIN_VIOL)
        {
            if (cp_add_lp_cut(cp, bac_env, cut))
            {
                rval = cp_add_lp_cut2data(cp, bac_env, cut, &data);
                check_rval(rval, " failed.", cleanup);

                if (data)
                {
                    tnadded++;
                    addpiece++;
                }

                if (data && data->nrows >= SOLVER_IP_BAC_STORE_BATCH)
                {
                    if (data->nrows > 0)
                    {
                        rval = lp_add_rows(lp, data);
                        check_rval(rval, " failed.", cleanup);
                        lp_free_data(&data);
                    }
                }

                if (addpiece > SOLVER_IP_BAC_CUT_BATCH)
                {
                    if (data->nrows > 0)
                    {
                        rval = lp_add_rows(lp, data);
                        check_rval(rval, " failed.", cleanup);
                        lp_free_data(&data);
                    }

                    rval = cp_update_lp(cp, bac_env, 0);
                    check_rval(rval, " failed.", cleanup);

                    if (lp->status == SOLVER_LP_INFEASIBLE)
                    {
                        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
                            printf("LP is really infeasible (process_cuts)\n");

                        while (bac_env->cut_queue)
                        {
                            cut                = bac_env->cut_queue;
                            bac_env->cut_queue = cut->next;
                            cp_free_cut(&cut);
                        }
                        goto cleanup;
                    }

                    addpiece = 0;
                }
            }
            else
            {
                cp_free_cut(&cut);
            }
        }
        else
        {
            cp_free_cut(&cut);
        }
    }

    if (addpiece > 0)
    {
        if (data)
        {
            rval = lp_add_rows(lp, data);
            check_rval(rval, " failed.", cleanup);
            lp_free_data(&data);
        }
        rval = cp_update_lp(cp, bac_env, 0);
        check_rval(rval, " failed.", cleanup);

        if (lp->status == SOLVER_LP_INFEASIBLE)
        {
            if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
                printf("LP is really infeasible (process_cuts2)\n");
            goto cleanup;
        }
    }

    if (oldval > lp->sol->val + SOLVER_ZEROPLUS)
    {
        rval = cp_update_lp(cp, bac_env, 1);
        check_rval(rval, "failed", cleanup);
    }

    rval = 0;

cleanup:

    *nadded = tnadded;
    return rval;
}

int
cp_add_lp_cut2data(cp_prob *cp, cp_exact_bac_env *bac_env, cp_cut *cut,
                   lp_data **data)
{
    int nzcnt;
    int i, j, end;
    int rval = 0;
    graph_arc **nzlist;
    graph_arc *arc;
    graph_vertex *v;

    lp_prob *lp         = bac_env->ip->lp;
    solver_graph *graph = lp->graph;

    if (!(*data))
    {
        *data = lp_create_data(0, SOLVER_IP_BAC_STORE_BATCH * 1.2, 10000);
        (*data)->nrows  = 0;
        (*data)->nzcnt  = 0;
        (*data)->beg[0] = 0;
    }

    if (cut->tcount + cut->hcount)
    {
        nzlist = malloc(graph->na * sizeof(graph_arc *));
        cp_get_cut_arcs(graph, cut, &nzcnt, nzlist);

        if (((*data)->nzcnt + nzcnt + graph->nv + 2) >= (*data)->nzspace)
            lp_realloc_data(*data, (*data)->nzcnt + nzcnt + graph->nv + 1, 2);

        for (i = 0, end = nzcnt, nzcnt = 0; i < end; i++)
        {
            arc = nzlist[i];
            if (arc->coef)
            {
                (*data)->val[(*data)->nzcnt + nzcnt] = arc->coef;
                (*data)->ind[(*data)->nzcnt + nzcnt] = graph->nv + arc->ind;
                arc->coef                            = 0;
                nzcnt++;
            }
        }

        for (i = 0; i < 2 * cut->tcount; i++)
        {
            (*data)->val[(*data)->nzcnt + nzcnt] = cut->vycoef;
            (*data)->ind[(*data)->nzcnt + nzcnt] = cut->verts[i];
            nzcnt++;
        }

        (*data)->rhs[(*data)->nrows]     = cut->rhs;
        (*data)->sense[(*data)->nrows++] = cut->sense;
        (*data)->nzcnt += nzcnt;
        free(nzlist);
    }
    else if (cut->logical)
    {
        if (((*data)->nzcnt + 2) >= (*data)->nzspace)
            lp_realloc_data(*data, (*data)->nzcnt + 2, 2);
        (*data)->val[(*data)->nzcnt]     = -1;
        (*data)->ind[(*data)->nzcnt]     = cut->logical->v;
        (*data)->val[(*data)->nzcnt + 1] = 1;
        arc = graph_find_arc_hash(graph->archash, cut->logical->arc[0],
                                  cut->logical->arc[1]);
        (*data)->ind[(*data)->nzcnt + 1] = graph->nv + arc->ind;

        (*data)->rhs[(*data)->nrows]     = 0.0;
        (*data)->sense[(*data)->nrows++] = 'L';
        (*data)->nzcnt += 2;
    }
    else if (cut->cover_edge)
    {
        if (((*data)->nzcnt + 2 * cut->cover_edge->na) >= (*data)->nzspace)
            lp_realloc_data(*data, (*data)->nzcnt + 2 * cut->cover_edge->na, 2);
        graph->marker++;
        for (i = 0; i < cut->cover_edge->na; i++)
        {
            (*data)->val[((*data)->nzcnt)] = 1;
            arc =
            graph_find_arc_hash(graph->archash, cut->cover_edge->arcs[2 * i],
                                cut->cover_edge->arcs[2 * i + 1]);
            (*data)->ind[((*data)->nzcnt)++] = graph->nv + arc->ind;

            if (cut->cover_edge->strong)
            {
                v = graph->v[cut->cover_edge->arcs[2 * i]];
                if (v->mark != graph->marker)
                {
                    (*data)->val[((*data)->nzcnt)]   = -1;
                    (*data)->ind[((*data)->nzcnt)++] = v->i;
                    v->mark                          = graph->marker;
                }

                v = graph->v[cut->cover_edge->arcs[2 * i + 1]];
                if (v->mark != graph->marker)
                {
                    (*data)->val[((*data)->nzcnt)]   = -1;
                    (*data)->ind[((*data)->nzcnt)++] = v->i;
                    v->mark                          = graph->marker;
                }
            }
        }

        (*data)->rhs[(*data)->nrows]     = cut->rhs;
        (*data)->sense[(*data)->nrows++] = 'L';
    }
    else if (cut->cover_vertex)
    {
        if (((*data)->nzcnt + clique_count(cut->cover_vertex->verts)) >=
            (*data)->nzspace)
            lp_realloc_data(
            *data, (*data)->nzcnt + clique_count(cut->cover_vertex->verts), 2);
        graph->marker++;
        FOREACH_NODE_IN_CLIQUE (i, cut->cover_vertex->verts, j)
            ;
        {
            (*data)->val[((*data)->nzcnt)]   = 1;
            (*data)->ind[((*data)->nzcnt)++] = i;
        }

        (*data)->rhs[(*data)->nrows]     = cut->rhs;
        (*data)->sense[(*data)->nrows++] = 'L';
    }
    else if (cut->path)
    {
        if (((*data)->nzcnt + 2 * cut->path->na + cut->path->fna) >=
            (*data)->nzspace)
            lp_realloc_data(
            *data, (*data)->nzcnt + 2 * cut->path->na + 2 * cut->path->fna, 2);
        for (i = 0; i < cut->path->na; i++)
        {
            (*data)->val[((*data)->nzcnt)] = 1;
            arc = graph_find_arc_hash(graph->archash, cut->path->arcs[2 * i],
                                      cut->path->arcs[2 * i + 1]);
            (*data)->ind[((*data)->nzcnt)++] = graph->nv + arc->ind;
        }
        for (i = 0; i < cut->path->fna; i++)
        {
            arc = graph_find_arc_hash(graph->archash, cut->path->farcs[2 * i],
                                      cut->path->farcs[2 * i + 1]);
            if (arc)
            {
                (*data)->val[((*data)->nzcnt)]   = -1;
                (*data)->ind[((*data)->nzcnt)++] = graph->nv + arc->ind;
            }
        }

        graph->marker++;
        for (i = 0; i < 2 * cut->path->na; i++)
        {
            if (graph->v[cut->path->arcs[i]]->mark == graph->marker)
            {
                (*data)->val[((*data)->nzcnt)]   = -1;
                (*data)->ind[((*data)->nzcnt)++] = cut->path->arcs[i];
            }
            graph->v[cut->path->arcs[i]]->mark = graph->marker;
        }

        (*data)->rhs[(*data)->nrows]     = cut->rhs;
        (*data)->sense[(*data)->nrows++] = 'L';
    }

    (*data)->beg[(*data)->nrows] = (*data)->nzcnt;

    return rval;
}

int
cp_update_lp(cp_prob *cp, cp_exact_bac_env *bac_env, int age)
{
    int rval;
    int ndeleted;
    double oldobjval, objval;

    lp_prob *lp = bac_env->ip->lp;

    oldobjval = lp->sol->val;

    rval = lp_opt_dual(lp);

    if (lp->status == SOLVER_LP_INFEASIBLE)
    {
        rval = cp_recover_infeas(cp, bac_env);
        check_rval(rval, " failed", cleanup);
        if (lp->status == SOLVER_LP_INFEASIBLE)
        {
            printf("Problem is really infeasible (cp_update_lp)\n");
            return 0;
        }
    }

    rval = cp_update_lp_sol(cp, bac_env);
    check_rval(rval, "failed", cleanup);

    objval = lp_get_objval(lp);

    if (age)
    {
        rval = stats_start(bac_env->stats->age_vars);
        check_rval(rval, "failed", cleanup);
        rval = cp_age_vars(cp, bac_env, &ndeleted);
        check_rval(rval, "failed", cleanup);
        if (ndeleted)
        {
            rval = lp_opt_dual(lp);
            if (lp->status == SOLVER_LP_INFEASIBLE)
                printf("AGE VARS LP INFEASIBLE\n");
            check_rval(rval, " failed", cleanup);
            rval = cp_update_lp_sol(cp, bac_env);
            check_rval(rval, "failed", cleanup);
        }
        rval = stats_stop(bac_env->stats->age_vars, ndeleted);
        check_rval(rval, "failed", cleanup);

        rval = stats_start(bac_env->stats->age_cuts);
        check_rval(rval, "failed", cleanup);
        check_rval(rval, "failed", cleanup);
        rval = cp_age_cuts(cp, bac_env, &ndeleted);
        if (ndeleted)
        {
            rval = lp_opt_dual(lp);
            if (lp->status == SOLVER_LP_INFEASIBLE)
                printf("AGE CUTS LP INFEASIBLE\n");
            check_rval(rval, " failed", cleanup);
            rval = cp_update_lp_sol(cp, bac_env);
            check_rval(rval, "failed", cleanup);
        }
        rval = stats_stop(bac_env->stats->age_cuts, ndeleted);
        check_rval(rval, "failed", cleanup);
    }

    return 0;

cleanup:

    return rval;
}

int
cp_update_lp_sol(cp_prob *cp, cp_exact_bac_env *bac_env)
{
    int rval = 0;
    graph_arc *arc, *lparc, *next;
    graph_vertex *v     = NULL;
    graph_arc **xsorted = NULL;
    double *x           = NULL;
    graph_vertex **heap = NULL;
    int i;
    int narcs;
    lp_prob *lp = bac_env->ip->lp;

    lp->sol->val      = lp_get_objval(lp);
    cp_cut_repo *cuts = bac_env->cuts;

    check_assert(lp->graph->nv + lp->graph->na == lp_get_ncols(lp), "",
                 CLEANUP);
    check_assert(lp->graph->nv + 1 + cuts->count == lp_get_nrows(lp), "",
                 CLEANUP);

    check_assert(lp->graph->nv + lp->graph->na == lp_get_ncols(lp), "",
                 CLEANUP);
    check_assert(cuts->count + 1 + lp->graph->nv == lp_get_nrows(lp), "",
                 CLEANUP);

    x = malloc((lp->graph->nv + lp->graph->na) * sizeof(double));
    check_null(x, "out of memory", CLEANUP);

    rval = lp_get_x(lp, x);
    check_rval(rval, " failed", CLEANUP);

    graph_free(&lp->sol->graph);
    lp->sol->graph       = graph_create();
    lp->sol->graph->data = cp->data;
    graph_add_vertices(lp->sol->graph, lp->graph->nv);
    for (i = 0; i < lp->graph->nv; i++)
    {
        lp->sol->graph->v[i]->fixed  = lp->graph->v[i]->fixed;
        lp->sol->graph->v[i]->branch = lp->graph->v[i]->branch;
        lp->sol->graph->v[i]->obj    = lp->graph->v[i]->obj;
    }
    lp->sol->graph->orig = lp->graph;

    lp->sol->graph->fixed = NULL;
    for (v = lp->graph->fixed; v; v = v->fixed_next)
    {
        if (lp->sol->graph->fixed)
            lp->sol->graph->fixed->fixed_prev = lp->sol->graph->v[v->i];
        lp->sol->graph->v[v->i]->fixed_next = lp->sol->graph->fixed;
        lp->sol->graph->fixed               = lp->sol->graph->v[v->i];
    }

    xsorted = malloc(lp->graph->na * sizeof(graph_arc *));
    check_null(xsorted, "out of memory\n", CLEANUP);

    lp->sol->integral = 1;
    for (i = 0, narcs = 0; i < lp->graph->na; i++)
    {
        lparc = lp->graph->arcs[i];
        if (x[lp->graph->nv + lparc->ind] > SOLVER_ZEROPLUS)
        {
            xsorted[narcs]      = lp->graph->arcs[i];
            xsorted[narcs++]->x = x[lp->sol->graph->nv + lparc->ind];
            if (x[lp->graph->nv + lparc->ind] <= SOLVER_ONEMINUS)
                lp->sol->integral = 0;
        }
    }

    for (i = 0; i < narcs; i++)
    {
        lparc = xsorted[i];
        arc   = graph_add_arc(lp->sol->graph, lparc->tail->i, lparc->head->i);
        arc->ind    = lparc->ind;
        arc->cost   = lparc->cost;
        arc->branch = lparc->branch;
        arc->x      = x[lp->sol->graph->nv + lparc->ind];
        arc->orig   = lparc;
        lparc->x    = 0.0;
        arc->tail->y += x[lp->sol->graph->nv + lparc->ind] / 2;
        arc->head->y += x[lp->sol->graph->nv + lparc->ind] / 2;
    }

    heap      = malloc(lp->sol->graph->n3v * sizeof(graph_vertex *));
    int nheap = 0;
    for (v = lp->sol->graph->tail; v; v = v->next)
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
                lp->sol->integral = 0;
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
                    graph_del_arc(lp->sol->graph, &arc);
                }
                lp->sol->integral = 0;
            }
            v->y = 0;
            x[i] = 0;
        }
    }

    if (lp->sol->integral)
    {
        int ncomp;
        int *compscount      = NULL;
        graph_vertex **comps = NULL;
        rval = graph_get_comps(lp->sol->graph, &ncomp, &compscount, &comps);
        check_rval(rval, "", CLEANUP);
        if (ncomp == 1)
            lp->sol->graph->connected = 1;
        free(compscount);
        free(comps);
    }

    if (lp->sol->x)
        free(lp->sol->x);

    lp->sol->x = x;

CLEANUP:
    if (rval && x)
        free(x);
    if (xsorted)
        free(xsorted);
    if (heap)
        free(heap);
    return rval;
}
