#include "cp/cp.h"
#include "cp/exact/bac/bac.h"

static int
test_cut(cp_prob *cp, cp_exact_bac_env *bac_env, solver_graph *graph,
         graph_clique_repo *repo, int cvcount, graph_vertex **cverts,
         double val);

int
cp_sep_sec_exact_hong(cp_prob *cp, cp_exact_bac_env *bac_env,
                      solver_graph *graph, graph_clique_repo *repo)
{
    int rval = 0;

    int reorder;

    double val;

    graph_vertex *squeue = NULL;
    graph_arc *f;
    graph_vertex **cverts = NULL;
    graph_vertex *other, *currvert, *nextvert;
    int cvcount;
    if (!graph->tail)
        goto CLEANUP;
    nextvert = graph->tail->next;
    while (nextvert)
    {
        currvert           = nextvert;
        nextvert->mark_aux = 0;

        if ((bac_env->param->sep_sec_exact & CP_SEC_HONG_ST) == CP_SEC_HONG_ST)
        {
            rval = graph_get_mincut_st(graph, graph->tail, currvert, &val,
                                       &cverts, &cvcount);
            check_rval(rval, " failed", CLEANUP);

            test_cut(cp, bac_env, graph, repo, cvcount, cverts, val);
            free(cverts);
        }

        if ((bac_env->param->sep_sec_exact & CP_SEC_HONG_TS) == CP_SEC_HONG_TS)
        {
            for (graph_vertex *v = graph->fixed; v; v = v->fixed_next)
            {
                if (currvert->i != v->i)
                {
                    check_assert(v->fixed, "", CLEANUP);
                    rval = graph_get_mincut_st(graph, currvert, v, &val,
                                               &cverts, &cvcount);
                    check_rval(rval, "solve_mincut failed", CLEANUP);

                    test_cut(cp, bac_env, graph, repo, cvcount, cverts, val);
                    free(cverts);
                }
            }
        }

        /* Shrink vertices if Rule S3 is considered */
        if (bac_env->param->srk_s3)
        {
            reorder = 0;

            f = graph_find_arc(graph, graph->tail, currvert);
            if (f && f->x > graph->tail->y)
                reorder = 1;

            graph_identify_vertices(graph, graph->tail, currvert);
            check_rval(rval, "identify_srkvertices failed", CLEANUP);

            if (graph->n3v > 1 && graph->na > 1)
            {
                if (reorder)
                    graph_reorder_vertices(graph);
                if (bac_env->param->srk_extra)
                {
                    squeue = NULL;
                    for (f = graph->tail->edge; f; f = outnext(f, graph->tail))
                    {
                        other        = otherend(f, graph->tail);
                        other->qnext = squeue;
                        squeue       = other;
                    }
                    graph->tail->qnext = squeue;
                    squeue             = graph->tail;

                    /****************************************************/

                    if (bac_env->param->srk_rule != CP_SRK_NONE &&
                        graph->na > 3)
                    {
                        rval = cp_shrink_exact_bac_graph(cp, bac_env, graph,
                                                         squeue, repo);
                        check_rval(rval, "shrink_cp_c1 failed", CLEANUP);
                    }
                }

                if (graph->n3v)
                {
                    nextvert = graph->tail->next;
                    while (nextvert && !(nextvert->mark_aux))
                        nextvert = nextvert->next;
                }
                else
                {
                    nextvert = NULL;
                }
            }
            else
            {
                nextvert = NULL;
            }
        }
        else
        {
            nextvert = graph->tail->next;
            while (nextvert && !(nextvert->mark_aux)) nextvert = nextvert->next;
        }
    }

CLEANUP:

    return rval;
}

static int
test_cut(cp_prob *cp, cp_exact_bac_env *bac_env, solver_graph *graph,
         graph_clique_repo *repo, int cvcount, graph_vertex **cverts,
         double val)
{
    int rval = 0;
    int i;
    graph_clique *clique = NULL;
    graph_vertex *v;
    graph_vertex **verts  = NULL;
    graph_vertex **tverts = NULL;
    int fixed_in          = 0;
    double obj_in         = 0;
    int nverts;
    double maxweight = 0.0;

    verts = malloc(graph->nv * sizeof(graph_vertex *));
    check_null(verts, "out of memory", CLEANUP);
    tverts = malloc(graph->n3v * sizeof(graph_vertex *));
    check_null(tverts, "out of memory", CLEANUP);

    if (val < 2.0 - SOLVER_ZEROPLUS)
    {
        for (i = 0, nverts = 0; i < cvcount; i++)
        {
            verts[nverts++] = cverts[i]->orig;
            obj_in += cverts[i]->orig->obj;
            if (cverts[i]->orig->fixed)
                fixed_in++;

            if (maxweight < cverts[i]->orig->y)
                maxweight = cverts[i]->orig->y;

            for (v = cverts[i]->members; v; v = v->members)
            {
                if (maxweight < v->orig->y)
                    maxweight = v->orig->y;
                verts[nverts++] = v->orig;
                obj_in += v->orig->obj;
                if (v->fixed)
                    fixed_in++;
            }
        }
        if (nverts <= graph->nv / 2)
            clique = clique_conv_vertices2clique(graph->orig, verts, nverts);
        else
            clique = clique_conv_vertices2coclique(graph->orig, verts, nverts);
        clique->val = val;
        clique_register_repo(graph->orig, repo, clique);
        clique_free(&clique);

        if (bac_env->param->sep_sec_cc_extra && cp->data->tot_obj_edge == 0.0 &&
            fixed_in)
        {
            if (bac_env->param->srk_rule != CP_SRK_NONE ||
                bac_env->param->srk_s3 != 0)
                cp_find_sec_cliques(graph, cverts, cvcount, val,
                                    cp->ip->lowerboundG, repo);
        }
        else if (bac_env->param->sep_sec_cc_extra &&
                 cp->data->tot_obj_edge == 0.0)
        {

            if (bac_env->param->srk_rule != CP_SRK_NONE ||
                bac_env->param->srk_s3 != 0)
            {
                graph->marker++;
                for (i = 0; i < cvcount; i++)
                {
                    check_assert(cverts[i]->mark < graph->marker, "", CLEANUP);
                    cverts[i]->mark = graph->marker;
                }

                int tcount = 0;
                for (v = graph->tail; v; v = v->next)
                {
                    if (v->mark != graph->marker && v->deg)
                        tverts[tcount++] = v;
                }

                cp_find_sec_cliques(graph, tverts, tcount, val,
                                    cp->ip->lowerboundG, repo);
            }
        }
    }

DONE:
    if (verts)
        free(verts);
    if (tverts)
        free(tverts);
CLEANUP:
    return rval;
}
