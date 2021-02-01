#include "cp/cp.h"
#include "ip/exact/bac/bac.h"
#include "cp/exact/bac/bac.h"

#define inhandle(arc, marker)                                                  \
    ((arc)->tail->mark == marker ? (arc)->tail : (arc)->head)
#define outhandle(arc, marker)                                                 \
    ((arc)->tail->mark != marker ? (arc)->tail : (arc)->head)

static int
grab_component(int *hcount, graph_vertex **handle, graph_vertex *v, double yval,
               int label),
grow_teeth(solver_graph *graph, int hcount, graph_vertex **handle, int *tcount,
           graph_arc **teeth, double yval, double *delta),
work_blossom(solver_graph *graph, int hcount, graph_vertex **handle, int tcount,
             graph_arc **teeth, double yval, int depth, cp_cut **cuts,
             int *cutcount, graph_clique_repo *repo),
add_blossom(solver_graph *graph, int hcount, graph_vertex **handle, int tcount,
            graph_arc **teeth, cp_cut **cuts, int *cutcount);

static int
sort_verts(const void *xx, const void *yy)
{
    graph_vertex *x = *(graph_vertex **)xx, *y = *(graph_vertex **)yy;

    if (x->y < y->y)
        return +1;
    if (x->y > y->y)
        return -1;

    return 0;
}

int
cp_sep_blossom_fast(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                    cp_cut **cuts)
{
    int rval               = 0;
    graph_vertex **vsort   = NULL;
    graph_vertex ***levels = NULL;
    graph_vertex **handle  = NULL;
    graph_arc **teeth      = NULL;
    graph_vertex *v;
    int i, j, k;
    int hcount, tcount, nlevels = 0;
    int *levelcnt = NULL;
    graph_clique *clique;
    double delta;
    graph_clique_repo *repo = NULL;
    cp_cut_list *cut_list   = NULL;
    lp_prob *lp             = bac_env->ip->lp;

    solver_graph *graph = lp->sol->graph;

    *cutcount = 0;
    *cuts     = NULL;

    vsort = malloc(graph->n3v * sizeof(graph_vertex *));
    check_null(vsort, "out of memory ", CLEANUP);

    repo = clique_create_repo(graph->nv);
    check_null(repo, "out of memory ", CLEANUP);
    cut_list = cp_create_cut_list(bac_env->param->sec_max_cut, 2.0);
    check_null(cut_list, "failed", CLEANUP);

    for (i = 0, k = 0; i < graph->nv; i++)
    {
        v = graph->v[i];
        if (v->deg)
            vsort[k++] = v;
        v->comp = 0;
    }

    assert(k == graph->n3v);

    qsort(vsort, graph->n3v, sizeof(graph_vertex *), sort_verts);

    levelcnt = malloc(graph->n3v * sizeof(int));
    check_null(levelcnt, "out of memory ", CLEANUP);
    memset(levelcnt, 0, graph->n3v * sizeof(int));

    levelcnt[0] = 1;
    for (i = 1, nlevels = 1;
         i < graph->n3v && vsort[i]->y > SOLVER_IP_BAC_MIN_VIOL; i++)
    {
        if (fabs(vsort[i - 1]->y - vsort[i]->y) > SOLVER_ZEROPLUS)
            nlevels++;
        levelcnt[nlevels - 1]++;
    }

    levels = malloc(k * sizeof(graph_vertex **));
    check_null(levels, "out of memory", CLEANUP);

    for (i = 0; i < nlevels; i++)
    {
        levels[i]   = malloc(levelcnt[i] * sizeof(graph_vertex *));
        levelcnt[i] = 0;
    }

    levelcnt[0]  = 1;
    levels[0][0] = vsort[0];
    for (i = 1, nlevels = 1;
         i < graph->n3v && vsort[i]->y > SOLVER_IP_BAC_MIN_VIOL; i++)
    {
        if (fabs(vsort[i - 1]->y - vsort[i]->y) > SOLVER_ZEROPLUS)
            nlevels++;
        levels[nlevels - 1][levelcnt[nlevels - 1]++] = vsort[i];
    }

    handle = malloc(graph->nv * sizeof(graph_vertex *));
    check_null(handle, "out of memory", CLEANUP);
    teeth = malloc(graph->na * sizeof(graph_arc *));
    check_null(teeth, "out of memory", CLEANUP);

    for (i = 0; i < nlevels; i++)
    {
        for (j = 0; j < graph->nv; j++) graph->v[j]->comp = 0;

        k = 0;
        for (j = 0; j < levelcnt[i]; j++)
        {
            v = levels[i][j];
            if (v->comp == 0)
            {
                hcount = 0;
                rval   = grab_component(&hcount, handle, v, v->y, ++k);
                check_rval(rval, "grab_component failed", CLEANUP);

                grow_teeth(graph, hcount, handle, &tcount, teeth, v->y, &delta);

                if (tcount > 2)
                {

                    rval = work_blossom(graph, hcount, handle, tcount, teeth,
                                        v->y, 0, cuts, cutcount, repo);
                    check_rval(rval, "failed", CLEANUP);
                }
                else if (tcount < 3 && delta < 2 && hcount >= 3)
                {
                    clique = clique_conv_vertices2clique(graph, handle, hcount);
                    clique->val = delta;
                    clique_register_repo(graph, repo, clique);
                    clique_free(&clique);
                    free(clique);
                }
            }
        }
    }

    cp_get_clique_repo_sec_cuts(cp, bac_env, graph, repo, cut_list);
    cp_get_cut_list(cp, cut_list, cutcount, cuts);

CLEANUP:

    if (handle)
        free(handle);
    if (teeth)
        free(teeth);

    if (vsort)
        free(vsort);
    for (i = 0; i < nlevels; i++) free(levels[i]);
    if (levels)
        free(levels);
    if (levelcnt)
        free(levelcnt);
    if (repo)
        clique_free_repo(&repo);
    if (cut_list)
        cp_free_cut_list(&cut_list);

    return rval;
}

static int
grab_component(int *hcount, graph_vertex **handle, graph_vertex *v, double yval,
               int label)
{
    graph_arc *arc;
    graph_vertex *other;

    v->comp             = label;
    handle[(*hcount)++] = v;

    for (arc = v->edge; arc; arc = outnext(arc, v))
    {
        if (arc->x > SOLVER_ZEROPLUS && arc->x < SOLVER_ONEMINUS * yval)
        {
            other = otherend(arc, v);
            if (other->comp == 0)
            {
                if (grab_component(hcount, handle, other, yval, label))
                    return 1;
            }
        }
    }
    return 0;
}

static int
grow_teeth(solver_graph *graph, int hcount, graph_vertex **handle, int *tcount,
           graph_arc **teeth, double yval, double *delta)
{
    int i;
    graph_vertex *v;
    graph_arc *arc;
    *tcount      = 0;
    int tmpcount = 0;

    *delta = 0;
    graph->marker++;
    for (i = 0; i < hcount; i++) handle[i]->mark = graph->marker;

    for (i = 0; i < hcount; i++)
    {
        v = handle[i];
        for (arc = v->edge; arc; arc = outnext(arc, v))
        {
            if (arc->tail->mark != arc->head->mark)
            {
                *delta += arc->x;
                if (arc->x >= yval * SOLVER_ONEMINUS)
                {
                    teeth[(tmpcount)++] = arc;
                }
            }
        }
    }

    *tcount = tmpcount;

CLEANUP:

    return 0;
}

static int
work_blossom(solver_graph *graph, int hcount, graph_vertex **handle, int tcount,
             graph_arc **teeth, double yval, int depth, cp_cut **cuts,
             int *cutcount, graph_clique_repo *repo)
{
    graph_vertex *in, *out;
    graph_arc **newteeth     = NULL;
    graph_vertex **newhandle = NULL;
    int i, newhcount, newtcount, rval = 0;
    int *del       = NULL;
    int *add       = NULL;
    int *hit       = NULL;
    int newblossom = 0;
    double delta;

    /* Clean up intersecting teeth */
    if (depth == 10)
        goto CLEANUP;

    hit = malloc(graph->nv * sizeof(int));
    check_null(hit, "out of memory ", CLEANUP);
    del = malloc(graph->nv * sizeof(int));
    check_null(del, "out of memory ", CLEANUP);
    add = malloc(graph->nv * sizeof(int));
    check_null(add, "out of memory ", CLEANUP);
    newhandle = malloc(graph->nv * sizeof(graph_vertex *));
    check_null(newhandle, "out of memory ", CLEANUP);
    newteeth = malloc(graph->na * sizeof(graph_arc *));
    check_null(newteeth, "out of memory ", CLEANUP);

    for (i = 0; i < graph->nv; i++)
    {
        add[i] = del[i] = hit[i] = 0;
        graph->v[i]->mark        = 0;
    }

    graph->marker++;
    for (i = 0; i < hcount; i++) handle[i]->mark = graph->marker;

    for (i = 0; i < tcount; i++)
    {
        in = inhandle(teeth[i], graph->marker);
        if (hit[in->i])
        {
            del[in->i] = 1;
            newblossom = 1;
        }
        else
            hit[in->i] = 1;
    }

    for (i = 0; i < tcount; i++)
    {
        out         = outhandle(teeth[i], graph->marker);
        hit[out->i] = 0;
    }
    for (i = 0; i < tcount; i++)
    {
        out = outhandle(teeth[i], graph->marker);
        in  = inhandle(teeth[i], graph->marker);
        if (!del[in->i])
        {
            if (hit[out->i])
            {
                add[out->i] = 1;
                newblossom  = 1;
            }
            else
                hit[out->i] = 1;
        }
    }

    for (i = 0; i < graph->nv; i++) hit[i] = 0;

    for (i = 0; i < hcount; i++) hit[handle[i]->i] = 1;

    for (i = 0; i < tcount; i++)
    {
        in = inhandle(teeth[i], graph->marker);
        if (del[in->i])
            hit[in->i] = 0;
        if (add[in->i])
            hit[in->i] = 1;
    }

    for (i = 0; i < tcount; i++)
    {
        out = outhandle(teeth[i], graph->marker);
        in  = inhandle(teeth[i], graph->marker);
        if (add[out->i])
            hit[out->i] = 1;
    }

    for (i = 0, newhcount = 0; i < graph->nv; i++)
    {
        if (hit[i])
            newhandle[newhcount++] = graph->v[i];
    }

    grow_teeth(graph, newhcount, newhandle, &newtcount, newteeth, yval, &delta);

    if (newtcount <= 2 && delta < 2 && newhcount >= 3)
    {
        graph_clique *clique;
        clique      = clique_conv_vertices2clique(graph, newhandle, newhcount);
        clique->val = delta;
        clique_register_repo(graph, repo, clique);
        clique_free(&clique);
        free(clique);
    }

    if (newblossom)
    {
        rval = work_blossom(graph, newhcount, newhandle, newtcount, newteeth,
                            yval, depth + 1, cuts, cutcount, repo);
        check_rval(rval, "work_blossom failed", CLEANUP);
    }
    else if (newhcount >= 3 && (newtcount % 2 == 1) && newtcount >= 3)
    {
        rval = add_blossom(graph, newhcount, newhandle, newtcount, newteeth,
                           cuts, cutcount);
        check_rval(rval, "add_blossom failed", CLEANUP);
    }

CLEANUP:

    if (hit)
        free(hit);
    if (del)
        free(del);
    if (add)
        free(add);
    if (newhandle)
        free(newhandle);
    if (newteeth)
        free(newteeth);

    return rval;
}

static int
sort_teeth(const void *xx, const void *yy)
/**************************************************************************/
{
    graph_arc *x = *(graph_arc **)xx, *y = *(graph_arc **)yy;

    if (x->tail->i < y->tail->i)
        return -1;
    if (x->tail->i > y->tail->i)
        return +1;

    return 0;
}

static int
add_blossom(solver_graph *graph, int hcount, graph_vertex **handle, int tcount,
            graph_arc **teeth, cp_cut **cuts, int *cutcount)
{
    int i, rval = 0;
    graph_vertex **verts;
    cp_cut *cut;

    verts = malloc(2 * sizeof(graph_vertex *));
    check_null(verts, "out of memory", CLEANUP);

    cut = cp_create_cut();
    check_null(cut, "out of memory", CLEANUP);

    cut->handles = malloc(sizeof(graph_clique *));
    check_null(cut->handles, "out of memory", CLEANUP);
    cut->teeth = malloc(tcount * sizeof(graph_clique *));
    check_null(cut->teeth, "out of memory", CLEANUP);
    cut->tcount = tcount;
    cut->hcount = 1;

    cut->handles[0] = clique_conv_vertices2clique(graph, handle, hcount);
    qsort(teeth, tcount, sizeof(graph_arc *), sort_teeth);

    for (i = 0; i < tcount; i++)
    {
        verts[0]      = teeth[i]->tail;
        verts[1]      = teeth[i]->head;
        cut->teeth[i] = clique_conv_vertices2clique(graph, verts, 2);
    }

    cut->verts = malloc(2 * tcount * sizeof(int));
    for (i = 0; i < tcount; i++)
    {
        cut->verts[2 * i]     = teeth[i]->tail->i;
        cut->verts[2 * i + 1] = teeth[i]->head->i;
    }
    cut->vycoef = -2;

    cut->rhs = 1 - tcount;

    cut->sense  = 'G';
    cut->branch = 0;

    rval = cp_build_cut_skeleton(cut, graph->nv);
    check_rval(rval, "", CLEANUP);

    cut->next = *cuts;
    *cuts     = cut;
    (*cutcount)++;

CLEANUP:

    if (rval)
    {
        if (cut)
        {
            cp_free_cut(&cut);
        }
    }
    free(verts);

    return rval;
}
