#include "op-solver.h"

int
cp_sep_sec_exact(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                 cp_cut **cuts)
{
    int rval               = 0;
    solver_graph *graph    = bac_env->ip->lp->sol->graph;
    solver_graph *srkgraph = NULL;

    graph_clique_repo *clique_repo = NULL;
    cp_cut_list *cut_list          = NULL;

    int stop = 0;
    int ncomp;
    int *compscount      = NULL;
    graph_vertex **comps = NULL;
    graph_clique *clique = NULL;

    *cutcount = 0;
    *cuts     = NULL;

    clique_repo = clique_create_repo(graph->nv);
    check_null(clique_repo, "failed", CLEANUP);
    cut_list = cp_create_cut_list(bac_env->param->sec_max_cut, 2.0);
    check_null(cut_list, "failed", CLEANUP);

    srkgraph = graph->shrunk = graph_create();
    rval                     = graph_copy(graph, srkgraph);
    check_rval(rval, "failed", CLEANUP);

    graph_reorder_vertices(srkgraph);

    if (bac_env->param->srk_rule != CP_SRK_NONE)
    {
        rval =
        cp_shrink_exact_bac_graph(cp, bac_env, srkgraph, NULL, clique_repo);
        check_rval(rval, "shrink_cp_c1 failed", CLEANUP);
    }

    rval = graph_get_comps(graph, &ncomp, &compscount, &comps);
    check_rval(rval, "", CLEANUP);

    for (int i = 0, k = 0; i < ncomp; k += compscount[i], i++)
    {
        int contains_zero = 0;
        double ymax       = 0.0;
        for (int j = 0; j < compscount[i]; j++)
        {
            if (comps[k + j]->fixed)
                contains_zero = 1;
            if (ymax < comps[k + j]->y)
                ymax = comps[k + j]->y;
        }
        if (!contains_zero && compscount[i] > 1)
        {

            if (compscount[i] <= graph->nv / 2)
                clique =
                clique_conv_vertices2clique(graph, comps + k, compscount[i]);
            else
                clique =
                clique_conv_vertices2coclique(graph, comps + k, compscount[i]);
            clique->val = 0.0;
            clique_register_repo(graph, clique_repo, clique);
            clique_free(&clique);

            for (int j = 0; j < compscount[i]; j++)
            {
                graph_arc *next;
                for (graph_arc *f = srkgraph->v[comps[k + j]->i]->edge; f;
                     f            = next)
                {
                    next = outnext(f, srkgraph->v[comps[k + j]->i]);
                    graph_del_arc(srkgraph, &f);
                }
            }

            if (1 || (compscount[i] > 2 && graph->nv - compscount[i] > 2))
            {
                if (2 * ymax > SOLVER_IP_BAC_MIN_VIOL)
                {
                    if (compscount[i] <= graph->nv / 2)
                        clique = clique_conv_vertices2clique(graph, comps + k,
                                                             compscount[i]);
                    else
                        clique = clique_conv_vertices2coclique(graph, comps + k,
                                                               compscount[i]);
                    clique->val = 0.0;
                    clique_register_repo(graph, clique_repo, clique);
                    clique_free(&clique);
                }
            }
            else if (compscount[i] == 2)
            {
                cp_cut *cut           = NULL;
                cp_cut_logic *logical = NULL;
                // Logical tail
                cut = cp_create_cut();
                check_null(cut, "", CLEANUP);
                cut->rhs   = 0.0;
                cut->sense = 'L';

                logical         = malloc(sizeof(cp_cut_logic));
                logical->arc[0] = comps[k]->i;
                logical->arc[1] = comps[k + 1]->i;
                logical->v      = comps[k]->i;

                cut->logical = logical;

                cp_insert_cut_list(cut_list, 2 * comps[k + 1]->y, cut, &stop);

                // Logical head
                cut = cp_create_cut();
                check_null(cut, "", CLEANUP);
                cut->rhs   = 0.0;
                cut->sense = 'L';

                logical         = malloc(sizeof(cp_cut_logic));
                logical->arc[0] = comps[k]->i;
                logical->arc[1] = comps[k + 1]->i;
                logical->v      = comps[k + 1]->i;

                cut->logical = logical;

                cp_insert_cut_list(cut_list, 2 * comps[k]->y, cut, &stop);
            }
        }
        else if (contains_zero && ncomp > 2)
        {
            if (compscount[i] <= graph->nv / 2)
                clique =
                clique_conv_vertices2clique(graph, comps + k, compscount[i]);
            else
                clique =
                clique_conv_vertices2coclique(graph, comps + k, compscount[i]);
            clique->val = 0.0;
            clique_register_repo(graph, clique_repo, clique);
            clique_free(&clique);
        }
    }

    if (srkgraph->tail)
    {

        if (((bac_env->param->sep_sec_exact & CP_SEC_HONG_ST) ==
             CP_SEC_HONG_ST) ||
            ((bac_env->param->sep_sec_exact & CP_SEC_HONG_TS) ==
             CP_SEC_HONG_TS))
        {
            rval = cp_sep_sec_exact_hong(cp, bac_env, srkgraph, clique_repo);
        }
        else if (bac_env->param->sep_sec_exact == CP_SEC_GOMORYHU)
        {
            rval =
            cp_sep_sec_exact_gomoryhu(cp, bac_env, srkgraph, clique_repo);
        }
    }

GET_CUTS:

    rval =
    cp_get_clique_repo_sec_cuts(cp, bac_env, graph, clique_repo, cut_list);
    check_rval(rval, "", CLEANUP);
    cp_get_cut_list(cp, cut_list, cutcount, cuts);

CLEANUP:
    if (clique_repo)
        clique_free_repo(&clique_repo);
    if (cut_list)
        cp_free_cut_list(&cut_list);
    if (srkgraph)
        graph_free(&srkgraph);
    if (comps)
        free(comps);
    if (compscount)
        free(compscount);
    graph->shrunk = NULL;
    return rval;
}

int
cp_sep_sec_comps(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                 cp_cut **cuts)
{
    int rval = 0;
    int i, j, k, ncomp;
    double score_sum_in, score_sum_out, score_sum_deg_out;
    double ymax_in, ymax_out;
    int fixed_in, fixed_out;
    int ninvert;
    graph_vertex **invert = NULL;
    graph_vertex **comps  = NULL;
    int *compscount       = NULL;
    graph_vertex *v;
    cp_cut *cut            = NULL;
    cp_cut_logic *logical  = NULL;
    cp_cut *logi_cuts      = NULL;
    cp_cut *first_logi_cut = NULL;
    int logi_cutcount      = 0;
    int stop;

    solver_graph *graph    = bac_env->ip->lp->sol->graph;
    graph->data            = cp->data;
    solver_graph *srkgraph = NULL;

    *cutcount = 0;
    *cuts     = NULL;

    graph_clique *clique           = NULL;
    graph_clique_repo *clique_repo = NULL;
    cp_cut_list *cut_list          = NULL;

    clique_repo = clique_create_repo(graph->nv);
    check_null(clique_repo, "failed", CLEANUP);
    cut_list = cp_create_cut_list(bac_env->param->sec_max_cut, 2.0);
    check_null(cut_list, "failed", CLEANUP);

    rval = graph_get_comps(graph, &ncomp, &compscount, &comps);
    check_rval(rval, "failed", CLEANUP);

    invert = malloc(graph->nv * sizeof(graph_vertex));

    for (i = 0, k = 0, stop = 0; i < ncomp; k += compscount[i], i++)
    {
        clique  = NULL;
        ymax_in = ymax_out = -1.0;
        fixed_in = fixed_out = 0;
        score_sum_in = score_sum_out = score_sum_deg_out = 0.0;
        ninvert                                          = 0;
        graph->marker++;

        for (j = 0; j < compscount[i]; j++)
        {

            check_assert(comps[k + j]->mark < graph->marker, "", CLEANUP);
            comps[k + j]->mark = graph->marker;
        }

        for (j = 0; j < graph->nv; j++)
        {
            v = graph->v[j];
            if (v->mark == graph->marker)
            {
                if (v->fixed)
                    fixed_in++;
                if (ymax_in < v->max->y)
                    ymax_in = v->max->y;
                score_sum_in += v->obj;
            }
            else
            {
                if (v->fixed)
                    fixed_out++;
                if (ymax_out < v->max->y)
                    ymax_out = v->max->y;
                score_sum_out += v->obj;
                if (v->deg)
                {
                    score_sum_deg_out += v->obj;
                    invert[ninvert++] = v;
                }
            }
        }

        int connectivity = 0;
        if (cp->data->tot_obj_edge == 0.0)
        {
            if (cp->ip->lowerboundG > 0.0 &&
                ((fixed_in && cp->ip->lowerboundG >= score_sum_in) ||
                 (fixed_out && cp->ip->lowerboundG >= score_sum_out)))
            {
                if (compscount[i] <= graph->nv / 2)
                    clique = clique_conv_vertices2clique(graph, comps + k,
                                                         compscount[i]);
                else
                    clique = clique_conv_vertices2coclique(graph, comps + k,
                                                           compscount[i]);
                check_null(clique, "failed", CLEANUP);
                clique_register_repo(graph, clique_repo, clique);
                clique_free(&clique);
                connectivity++;
            }

            if (fixed_out && cp->ip->lowerboundG > 0.0 &&
                cp->ip->lowerboundG >= score_sum_deg_out)
            {
                clique = clique_conv_vertices2clique(graph, invert, ninvert);
                check_null(clique, "failed", CLEANUP);
                clique_register_repo(graph, clique_repo, clique);
                clique_free(&clique);
            }
        }

        if (!connectivity && !fixed_in)
        {
            if (1 || (compscount[i] > 2 && graph->nv - compscount[i] > 2))
            {
                if (2 * ymax_in + 2 * ymax_out > 2.0 + SOLVER_IP_BAC_MIN_VIOL)
                {
                    if (1 || compscount[i] <= graph->nv / 2)
                        clique = clique_conv_vertices2clique(graph, comps + k,
                                                             compscount[i]);
                    else
                        clique = clique_conv_vertices2coclique(graph, comps + k,
                                                               compscount[i]);
                    check_null(clique, "failed", CLEANUP);
                    clique_register_repo(graph, clique_repo, clique);
                    clique_free(&clique);
                }
            }
            else if (compscount[i] == 2)
            {
                // Logical tail
                cut = cp_create_cut();
                check_null(cut, "", CLEANUP);
                cut->rhs   = 0.0;
                cut->sense = 'L';

                logical         = malloc(sizeof(cp_cut_logic));
                logical->arc[0] = comps[k]->i;
                logical->arc[1] = comps[k + 1]->i;
                logical->v      = comps[k]->i;

                cut->logical = logical;

                cut->next = logi_cuts;
                logi_cuts = cut;
                logi_cutcount++;
                if (!first_logi_cut)
                    first_logi_cut = cut;

                // Logical head
                cut = cp_create_cut();
                check_null(cut, "", CLEANUP);
                cut->rhs   = 0.0;
                cut->sense = 'L';

                logical         = malloc(sizeof(cp_cut_logic));
                logical->arc[0] = comps[k]->i;
                logical->arc[1] = comps[k + 1]->i;
                logical->v      = comps[k + 1]->i;

                cut->logical = logical;

                cut->next = logi_cuts;
                logi_cuts = cut;
                logi_cutcount++;
            }
        }
    }

    rval =
    cp_get_clique_repo_sec_cuts(cp, bac_env, graph, clique_repo, cut_list);
    check_rval(rval, "", CLEANUP);
    cp_get_cut_list(cp, cut_list, cutcount, cuts);

    if (logi_cutcount)
    {
        first_logi_cut->next = *cuts;
        *cuts                = logi_cuts;
        cutcount += logi_cutcount;
    }

CLEANUP:
    if (invert)
        free(invert);
    if (comps)
        free(comps);
    if (compscount)
        free(compscount);
    if (clique_repo)
        clique_free_repo(&clique_repo);
    if (cut_list)
        cp_free_cut_list(&cut_list);
    graph->shrunk = NULL;

    return rval;
}

static void
choose(graph_vertex **dest, int k, graph_vertex **src, int n)
{
    int i, j = 0;

    assert(k <= n);

    for (i = 0; i < n && j < k; i++)
    {
        if ((rand() % (n - i)) < k - j)
        {
            dest[j] = src[i];
            j++;
        }
    }
}

static int
sort_vert(const void *xx, const void *yy)
/**************************************************************************/
{
    graph_vertex *x = *(graph_vertex **)xx, *y = *(graph_vertex **)yy;

    if (x->fixed)
        return -1;
    if (y->fixed)
        return +1;

    if (x->y > y->y)
        return +1;
    if (x->y < y->y)
        return -1;

    if (x->obj < y->obj)
        return -1;
    if (x->obj > y->obj)
        return +1;

    return 0;
}

int
cp_find_sec_cliques(solver_graph *srkgraph, graph_vertex **sverts, int scount,
                    double cutval, double lowerbound, graph_clique_repo *repo)
{
    int rval = 0;
    graph_vertex *v;
    graph_vertex **verts = NULL;
    graph_vertex **orig  = NULL;
    int contains_zero, orig_zero;
    int i, j;
    double tot_obj, tmpobj, orig_obj;
    double delta;
    int cnt, origcnt;
    int greedycount;
    graph_clique *clique;

    verts = malloc(srkgraph->nv * sizeof(graph_vertex *));
    check_null(verts, "out of memory\n", CLEANUP);

    qsort(sverts, scount, sizeof(graph_vertex *), sort_vert);

    // STRATEGY 1
#if 1
    if (sverts[0]->y < 1 - SOLVER_ZEROPLUS && sverts[0]->obj <= lowerbound)
    {
        cnt           = 0;
        contains_zero = 0;
        verts[cnt++]  = sverts[0]->orig;
        tmpobj        = sverts[0]->orig->obj;
        if (sverts[0]->orig->i == 0)
            contains_zero = 1;
        for (v = sverts[0]->members; v; v = v->members)
        {
            if (v->orig->i == 0)
                contains_zero = 1;
            tmpobj += v->orig->obj;
            verts[cnt++] = v->orig;
        }
        assert(tmpobj <= lowerbound);
        assert(contains_zero);

        if (cnt <= srkgraph->nv / 2)
            clique = clique_conv_vertices2clique(srkgraph->orig, verts, cnt);
        else
            clique = clique_conv_vertices2coclique(srkgraph->orig, verts, cnt);
        clique->val = 2 * sverts[0]->y;
        clique_register_repo(srkgraph->orig, repo, clique);
        clique_free(&clique);
    }
#endif

    for (i = 0, tot_obj = 0.0; i < scount; i++) tot_obj += sverts[i]->obj;

#if 1
    // STRATEGY 2
    orig = malloc(srkgraph->nv * sizeof(graph_vertex *));
    check_null(orig, "out of memory\n", CLEANUP);
    origcnt         = 0;
    cnt             = 0;
    orig[origcnt++] = sverts[0]->orig;
    verts[cnt++]    = sverts[0]->orig;

    contains_zero  = 0;
    double testobj = sverts[0]->obj;
    tmpobj         = sverts[0]->orig->obj;
    if (sverts[0]->orig->i == 0)
        contains_zero = 1;
    for (v = sverts[0]->members; v; v = v->members)
    {
        orig[origcnt++] = v->orig;
        verts[cnt++]    = v->orig;
        tmpobj += v->orig->obj;
        if (v->orig->i == 0)
            contains_zero = 1;
    }
    assert(contains_zero);
    assert(testobj == tmpobj);

    for (i = 1;
         i < scount && cutval + 2 * sverts[i]->y <= 2.0 - SOLVER_ZEROPLUS; i++)
    {
        orig_zero = contains_zero;
        orig_obj  = tmpobj;
        if (tot_obj - sverts[i]->obj <= lowerbound)
        {
            if (i < scount - 1)
            {
                for (j = i + 1; j < scount; j++)
                {
                    verts[cnt++] = sverts[j]->orig;
                    tmpobj += sverts[j]->orig->obj;
                    if (sverts[j]->orig->i == 0)
                        contains_zero = 1;
                    for (v = sverts[j]->members; v; v = v->members)
                    {
                        verts[cnt++] = v->orig;
                        tmpobj += v->orig->obj;
                        if (v->orig->i == 0)
                            contains_zero = 1;
                    }
                }
            }

            assert(tmpobj <= lowerbound);
            assert(contains_zero);

            if (cnt <= srkgraph->nv / 2)
                clique =
                clique_conv_vertices2clique(srkgraph->orig, verts, cnt);
            else
                clique =
                clique_conv_vertices2coclique(srkgraph->orig, verts, cnt);
            clique->val = cutval + 2 * sverts[i]->y;
            clique_register_repo(srkgraph->orig, repo, clique);
            clique_free(&clique);
        }

        if (i < scount - 1)
        {
            for (j = 0, cnt = 0; j < origcnt; j++) verts[cnt++] = orig[j];
            tmpobj          = orig_obj;
            contains_zero   = orig_zero;
            orig[origcnt++] = sverts[i]->orig;
            verts[cnt++]    = sverts[i]->orig;
            tmpobj += sverts[i]->orig->obj;
            if (sverts[i]->orig->i == 0)
                contains_zero = 1;
            for (v = sverts[i]->members; v; v = v->members)
            {
                orig[origcnt++] = v->orig;
                verts[cnt++]    = v->orig;
                tmpobj += v->orig->obj;
                if (v->orig->i == 0)
                    contains_zero = 1;
            }
        }
    }
#endif

#if 1
    // STRATEGY 3
    for (i = 1, greedycount = 0, delta = cutval, tmpobj = tot_obj;
         i < scount && delta + 2 * sverts[i]->y < 2.0 - SOLVER_ZEROPLUS; i++)
    {
        delta += 2 * sverts[i]->y;
        tmpobj -= sverts[i]->obj;
        greedycount++;
    }

    if (tmpobj <= lowerbound)
    {
        contains_zero = 0;
        cnt           = 0;
        verts[cnt++]  = sverts[0]->orig;
        tmpobj        = sverts[0]->orig->obj;
        if (sverts[0]->orig->i == 0)
            contains_zero = 1;
        for (v = sverts[0]->members; v; v = v->members)
        {
            tmpobj += v->orig->obj;
            verts[cnt++] = v->orig;
            if (v->orig->i == 0)
                contains_zero = 1;
        }
        for (i = scount - 1; i > greedycount; i--)
        {
            verts[cnt++] = sverts[i]->orig;
            tmpobj += sverts[i]->orig->obj;
            if (sverts[i]->orig->i == 0)
                contains_zero = 1;
            for (v = sverts[i]->members; v; v = v->members)
            {
                tmpobj += v->orig->obj;
                verts[cnt++] = v->orig;
                if (v->orig->i == 0)
                    contains_zero = 1;
            }
        }

        assert(tmpobj <= lowerbound);
        assert(contains_zero);

        if (cnt <= srkgraph->nv / 2)
            clique = clique_conv_vertices2clique(srkgraph->orig, verts, cnt);
        else
            clique = clique_conv_vertices2coclique(srkgraph->orig, verts, cnt);
        clique->val = delta;
        clique_register_repo(srkgraph->orig, repo, clique);
        clique_free(&clique);
    }
#endif

    rval = 0;

CLEANUP:

    if (orig)
        free(orig);
    if (verts)
        free(verts);

    return rval;
}

static int
sort_vert_y(const void *xx, const void *yy)
/**************************************************************************/
{
    graph_vertex *x = *(graph_vertex **)xx, *y = *(graph_vertex **)yy;

    if (x->fixed)
        return -1;
    if (y->fixed)
        return +1;

    if (x->y > y->y)
        return -1;
    if (x->y < y->y)
        return +1;

    if (x->obj > y->obj)
        return -1;
    if (x->obj < y->obj)
        return +1;

    return 0;
}

//#define SUPPOT_COCLIQUE

int
cp_get_clique_repo_sec_cuts(cp_prob *cp, cp_exact_bac_env *bac_env,
                            solver_graph *graph, graph_clique_repo *repo,
                            cp_cut_list *cut_list)
{
    int rval = 0;
    graph_clique *clique, *clique_orig;
    cp_cut *cut = NULL;

    int i, j, k;
    int nin = 0, nout = 0;
    int nselin = 0, nselout = 0;
    int stop, tmpstop;
    double cut_val;
    graph_clique *cliques[2]; //= {NULL};

    double ymax_in;
    double ymax_out;
    int fixed_in  = 0;
    int fixed_out = 0;
    int min_depth_in;
    int min_depth_out;
    double score_sum_in, score_sum_out;
    cp_cut_list *tmp_list = NULL;

    graph_vertex **max_in  = NULL;
    graph_vertex **max_out = NULL;
    graph_vertex **sel_in  = NULL;
    graph_vertex **sel_out = NULL;

    tmp_list = cp_create_cut_list(bac_env->param->sec_max_cut_x_clique, 2.0);
    check_null(tmp_list, "failed", CLEANUP);

#ifdef SUPPOT_COCLIQUE
    graph_vertex *v;
    graph_vertex **tverts = NULL;
    tverts                = malloc(graph->n3v * sizeof(graph_vertex *));
    graph_clique *clique_co;
    cp_cut_list *tmp_list2 = NULL;
    tmp_list2 = cp_create_cut_list(bac_env->param->sec_max_cut_x_clique, 2.0);
    check_null(tmp_list2, "failed", CLEANUP);
#endif

    max_in = malloc(graph->nv * sizeof(graph_vertex *));
    check_null(max_in, "out of memory in repo_clique_get_cuts", CLEANUP);
    max_out = malloc(graph->nv * sizeof(graph_vertex *));
    check_null(max_out, "out of memory in repo_clique_get_cuts", CLEANUP);

    sel_in = malloc(graph->n3v * sizeof(graph_vertex *));
    check_null(sel_in, "out of memory in repo_clique_get_cuts", CLEANUP);
    sel_out = malloc(graph->n3v * sizeof(graph_vertex *));
    check_null(sel_out, "out of memory in repo_clique_get_cuts", CLEANUP);

    for (i = 0, stop = 0; i < repo->size && !stop; i++)
    {
        clique_orig = repo->cliques[i];
        cliques[0]  = clique_orig;
#ifdef SUPPOT_COCLIQUE
        graph->marker++;
        FOREACH_NODE_IN_CLIQUE (k, clique_orig, j)
        {
            check_assert(graph->v[k]->mark < graph->marker, "", CLEANUP);
            graph->v[k]->mark = graph->marker;
        }

        int tcount = 0;
        for (v = graph->tail; v; v = v->next)
        {
            if (v->mark != graph->marker)
                tverts[tcount++] = v;
        }
        clique_co  = clique_conv_vertices2clique(graph, tverts, tcount);
        cliques[1] = clique_co;
#endif

#ifdef SUPPOT_COCLIQUE
        for (int l = 0; l < 2; l++)
#else
        for (int l = 0; l < 1; l++)
#endif
        {
            clique = cliques[l];

            ymax_in = ymax_out = -1.0;
            nin = nout = 0;
            fixed_in = fixed_out = 0;
            score_sum_in = score_sum_out = 0.0;
            min_depth_in = min_depth_out = SOLVER_MAXINT;

            graph->marker++;
            FOREACH_NODE_IN_CLIQUE (k, clique, j)
            {
                check_assert(graph->v[k]->mark < graph->marker, "", CLEANUP);
                graph->v[k]->mark = graph->marker;
            }

            for (k = 0; k < graph->nv; k++)
            {
                if (graph->v[k]->mark == graph->marker)
                {
                    score_sum_in += graph->v[k]->obj;
                    if (graph->v[k]->fixed)
                    {
                        if (!fixed_in)
                            nin = 0;
                        ymax_in       = graph->v[k]->y;
                        max_in[nin++] = graph->v[k];
                        fixed_in      = 1;
                        if (graph->v[k]->branch < min_depth_in)
                            min_depth_in = graph->v[k]->branch;
                    }
                    else if (!fixed_in)
                    {
                        if (ymax_in < graph->v[k]->y)
                        {
                            nin           = 0;
                            ymax_in       = graph->v[k]->y;
                            max_in[nin++] = graph->v[k];
                        }
                        else if (ymax_in < graph->v[k]->y + SOLVER_ZEROPLUS)
                        {
                            ymax_in       = graph->v[k]->y;
                            max_in[nin++] = graph->v[k];
                        }
                    }
                }
                else
                {
                    score_sum_out += graph->v[k]->obj;
                    if (graph->v[k]->fixed)
                    {
                        if (!fixed_out)
                            nout = 0;
                        ymax_out        = graph->v[k]->y;
                        max_out[nout++] = graph->v[k];
                        fixed_out       = 1;
                        if (graph->v[k]->branch < min_depth_out)
                            min_depth_out = graph->v[k]->branch;
                    }
                    else if (!fixed_out)
                    {
                        if (ymax_out < graph->v[k]->y)
                        {
                            nout            = 0;
                            ymax_out        = graph->v[k]->y;
                            max_out[nout++] = graph->v[k];
                        }
                        else if (ymax_out < graph->v[k]->y + SOLVER_ZEROPLUS)
                        {
                            ymax_out        = graph->v[k]->y;
                            max_out[nout++] = graph->v[k];
                        }
                    }
                }
            }

            if ((cp->data->tot_obj_edge == 0.0) &&
                ((fixed_in && score_sum_in <= cp->ip->lowerboundG) ||
                 (fixed_out && score_sum_out <= cp->ip->lowerboundG)))
            {
                cut = cp_create_cut();
                check_null(cut, "out of memory in conv_subtour0cuts", CLEANUP);

                cut->hcount = 1;
                cut->tcount = 0;

                cut->handles = malloc(sizeof(graph_clique *));
                check_null(cut->handles, "out of memory", CLEANUP);
                cut->handles[0] = clique_create();
                check_null(cut->handles[0], "out of memory", CLEANUP);
                clique_copy(clique, cut->handles[0]);

                cut->rhs   = 2.0;
                cut->sense = 'G';

                cp_insert_cut_list(cut_list, 2.0 - clique->val, cut, &tmpstop);
            }
            else
            {

                if (bac_env->param->sec_max_viol)
                {

                    nselin = bac_env->param->sec_max_vin < nin
                             ? bac_env->param->sec_max_vin
                             : nin;
                    nselout = bac_env->param->sec_max_vout < nout
                              ? bac_env->param->sec_max_vout
                              : nout;

                    choose(sel_out, nselout, max_out, nout);
                }
                else
                {
#if 1
                    for (k = 0, nin = 0, nout = 0;
                         k < graph->nv && fixed_in + fixed_out < 2; k++)
                    {
                        if (!fixed_in && graph->v[k]->mark == graph->marker)
                        {
                            if (clique->val - 2 * ymax_out + 2 +
                                SOLVER_IP_BAC_MIN_VIOL <
                                2 * graph->v[k]->y)
                            {
                                sel_in[nin++] = graph->v[k];
                            }
                        }
                        else if (!fixed_out &&
                                 graph->v[k]->mark != graph->marker)
                        {
                            if (clique->val - 2 * ymax_in + 2 +
                                SOLVER_IP_BAC_MIN_VIOL <
                                2 * graph->v[k]->y)
                            {
                                sel_out[nout++] = graph->v[k];
                            }
                        }
                    }
#endif

                    if (fixed_in)
                    {
                        sel_in[0] = max_in[0];
                        nin       = 1;
                    }
                    else
                        qsort(sel_in, nin, sizeof(graph_vertex *), sort_vert_y);
                    if (fixed_out)
                    {
                        sel_out[0] = max_out[0];
                        nout       = 1;
                    }
                    else
                        qsort(sel_out, nout, sizeof(graph_vertex *),
                              sort_vert_y);
                }

                // If a node d in the outside of the clique is fixed, i.e
                // d->y=1, then this can be simplified by max_out[0]=d, nout =
                // 1;

                for (k = 0; k < nout && k < bac_env->param->sec_max_vout; k++)
                {
                    if (bac_env->param->sec_max_viol)
                        choose(sel_in, nselin, max_in, nin);
                    for (j = 0; j < nin && j < bac_env->param->sec_max_vin; j++)
                    {
                        cut_val =
                        2 * sel_in[j]->y + 2 * sel_out[k]->y - 2 - clique->val;
                        if (cut_val > SOLVER_IP_BAC_MIN_VIOL)
                        {
                            // A
                            cut = cp_create_cut();
                            check_null(cut, "out of memory", CLEANUP);

                            cut->hcount = 0;
                            cut->tcount = 1;

                            cut->teeth = malloc(sizeof(graph_clique *));
                            check_null(cut->teeth, "out of memory", CLEANUP);
                            cut->teeth[0] = clique_create();
                            check_null(cut->teeth[0], "out of memory", CLEANUP);
                            clique_copy(clique, cut->teeth[0]);
                            cut->verts    = malloc(2 * sizeof(int));
                            cut->verts[0] = sel_in[j]->i;
                            cut->verts[1] = sel_out[k]->i;
                            cut->vycoef   = -2.0;

                            cut->rhs   = -2.0;
                            cut->sense = 'G';

                            check_assert(sel_in[j]->mark == graph->marker, "",
                                         CLEANUP);
                            check_assert(sel_out[k]->mark != graph->marker, "",
                                         CLEANUP);

                            cp_insert_cut_list(tmp_list, cut_val, cut,
                                               &tmpstop);
                        }
                    }
                }
            }
            cp_append_cut_list(cut_list, tmp_list, &stop);
            cp_erase_cut_list(tmp_list);
        }
    }

CLEANUP:

#ifdef SUPPOT_COCLIQUE
    if (tverts)
        free(tverts);
#endif
    if (max_in)
        free(max_in);
    if (max_out)
        free(max_out);
    if (sel_in)
        free(sel_in);
    if (sel_out)
        free(sel_out);
    if (tmp_list)
        cp_free_cut_list(&tmp_list);

    return rval;
}
