#include "cp/cp.h"

#define inhandle(arc, marker)                                                  \
    ((arc)->tail->mark == marker ? (arc)->tail : (arc)->head)
#define outhandle(arc, marker)                                                 \
    ((arc)->tail->mark != marker ? (arc)->tail : (arc)->head)

#define BLOTOLERANCE .01

#define GHFAST_EPSILON 0.3

static int
grab_component(int *hcount, graph_vertex **handle, graph_vertex *v, double yval,
               int label),
sort_teeth(const void *xx, const void *yy),
grow_ghteeth(solver_graph *graph, graph_clique_repo *repo, int hcount,
             graph_vertex **handle, int *tcount, graph_arc **teeth, double yval,
             int *cutcount, cp_cut **cuts),
work_blossom(solver_graph *graph, graph_clique_repo *repo, int hcount,
             graph_vertex **handle, int tcount, graph_arc **teeth, double yval,
             int depth, int *cutcount, cp_cut **cuts),
add_blossom(solver_graph *graph, int hcount, graph_vertex **handle, int tcount,
            graph_arc **teeth, cp_cut **cuts, int *cutcount);

static int
sort_verts(const void *xx, const void *yy)
/**************************************************************************/
{
    graph_vertex *x = *(graph_vertex **)xx, *y = *(graph_vertex **)yy;

    if (x->y < y->y)
        return +1;
    if (x->y > y->y)
        return -1;

    return 0;
}

int
cp_sep_blossom_ghfast(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                      cp_cut **cuts)
{
    int rval               = 0;
    graph_vertex **handle  = NULL;
    graph_arc **teeth      = NULL;
    graph_vertex **vsort   = NULL;
    graph_vertex ***levels = NULL;
    cp_cut_list *cut_list  = NULL;
    int *levelcnt          = NULL;
    int nlevels            = 0;
    graph_vertex *v;
    int i, j, k;
    int hcount, tcount;
    graph_clique_repo *repo = NULL;

    lp_prob *lp         = bac_env->ip->lp;
    solver_graph *graph = lp->sol->graph;

    *cutcount = 0;
    *cuts     = NULL;

    repo = clique_create_repo(graph->nv);
    check_null(repo, "out of memory", CLEANUP);
    cut_list = cp_create_cut_list(bac_env->param->sec_max_cut, 2.0);
    check_null(cut_list, "failed", CLEANUP);

    vsort = malloc(graph->n3v * sizeof(graph_vertex *));
    check_null(vsort, "out of memory", CLEANUP);

    for (i = 0, k = 0; i < graph->nv; i++)
    {
        v = graph->v[i];
        if (v->deg)
        {
            vsort[k++] = v;
        }
        v->comp = 0;
    }
    assert(k == graph->n3v);

    qsort(vsort, graph->n3v, sizeof(graph_vertex *), sort_verts);

    levelcnt = malloc(graph->n3v * sizeof(int));
    check_null(levelcnt, "out of memory", CLEANUP);
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

    for (i = 0; i < graph->nv; i++) graph->v[i]->comp = 0;

    handle = malloc(graph->nv * sizeof(graph_vertex *));
    teeth  = malloc(graph->na * sizeof(graph_arc *));

    for (i = 0; i < nlevels; i++)
    {
        for (j = 0; j < graph->nv; j++)
        {
            graph->v[j]->comp = 0;
            graph->v[j]->mark = 0;
        }

        k = 0;
        for (j = 0; j < levelcnt[i]; j++)
        {
            v = levels[i][j];
            if (v->comp == 0)
            {
                hcount = 0;
                rval   = grab_component(&hcount, handle, v, v->y, ++k);
                check_rval(rval, "grab_component failed", CLEANUP);

                grow_ghteeth(graph, repo, hcount, handle, &tcount, teeth, v->y,
                             cutcount, cuts);
            }
        }
    }

    cp_get_clique_repo_sec_cuts(cp, bac_env, graph, repo, cut_list);
    cp_get_cut_list(cp, cut_list, cutcount, cuts);

CLEANUP:

    for (i = 0; i < nlevels; i++) free(levels[i]);
    free(levels);
    free(levelcnt);
    free(vsort);

    free(handle);
    free(teeth);
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
    graph_arc *arc, *next;
    graph_vertex *other;

    v->comp             = label;
    handle[(*hcount)++] = v;

    for (arc = v->edge; arc; arc = next)
    {
        next = outnext(arc, v);
        if (arc->x >= GHFAST_EPSILON && arc->x <= yval * (1.0 - GHFAST_EPSILON))
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
sort_teeth(const void *xx, const void *yy)
/**************************************************************************/
{
    graph_arc *x = *(graph_arc **)xx, *y = *(graph_arc **)yy;

    if (x->x < y->x)
        return -1;
    if (x->x > y->x)
        return +1;

    return 0;
}

static int
grow_ghteeth(solver_graph *graph, graph_clique_repo *repo, int hcount,
             graph_vertex **handle, int *tcount, graph_arc **teeth, double yval,
             int *cutcount, cp_cut **cuts)
{
    int i, rval = 0;
    graph_vertex *v, *other;
    graph_arc *edge, *next;
    graph_arc *emax = NULL;
    *tcount         = 0;
    double z = 0.0, xemax = 0.0;
    double yhandle = 0.0;
    double delta;

    graph->marker++;

    for (i = 0; i < hcount; i++)
    {
        yhandle += handle[i]->y;
        handle[i]->mark = graph->marker;
    }

    delta = 0.0;
    for (i = 0; i < hcount; i++)
    {
        v = handle[i];
        for (edge = v->edge; edge; edge = next)
        {
            next = outnext(edge, v);
            if (edge->x > SOLVER_ZEROPLUS)
            {
                other = otherend(edge, v);
                if (other->mark != v->mark)
                {
                    delta += edge->x;
                    if (edge->x > (1.0 - GHFAST_EPSILON) * yval)
                    {
                        teeth[(*tcount)++] = edge;
                    }
                    else if (edge->x < GHFAST_EPSILON)
                    {
                        if (edge->x > xemax)
                        {
                            xemax = edge->x;
                            emax  = edge;
                        }
                    }
                }
                else
                {
                    z += edge->x;
                }
            }
        }
    }
    z *= 0.5;

    if (*tcount == 0)
    {
        goto CLEANUP;
    }
    else if (*tcount < 3 && hcount > 1 && delta < 2.)
    {

        graph_clique *clique =
        clique_conv_vertices2clique(graph, handle, hcount);
        clique->val = delta;
        clique_register_repo(graph, repo, clique);
        clique_free(&clique);
        free(clique);

        if (hcount < 3)
            goto CLEANUP;
    }
    else if (*tcount % 2 == 0)
    {
        if (emax)
            teeth[(*tcount)++] = emax;
    }

    qsort(teeth, *tcount, sizeof(graph_arc *), sort_teeth);

    if (*tcount % 2 == 0)
        (*tcount)--;

    i = 0;
    z += teeth[i++]->x;
    if (z > yhandle + (double)((i - 1) / 2) + BLOTOLERANCE)
    {
        rval = work_blossom(graph, repo, hcount, handle, i, teeth, yval, 0,
                            cutcount, cuts);
        check_rval(rval, "work_blossom failed", CLEANUP);
    }
    else
    {
        while (i < *tcount)
        {
            z += teeth[i++]->x;
            z += teeth[i++]->x;

            if (z > yhandle + (double)((i - 1) / 2) + BLOTOLERANCE)
            {
                rval = work_blossom(graph, repo, hcount, handle, i, teeth, yval,
                                    0, cutcount, cuts);
                check_rval(rval, "work_blossom failed", CLEANUP);
                break;
            }
        }
    }

CLEANUP:
    return rval;
}

static int
work_blossom(solver_graph *graph, graph_clique_repo *repo, int hcount,
             graph_vertex **handle, int tcount, graph_arc **teeth, double yval,
             int depth, int *cutcount, cp_cut **cuts)
{
    graph_vertex *in, *out;
    graph_arc **newteeth     = NULL;
    graph_vertex **newhandle = NULL;
    int i, newhcount, newtcount, rval = 0;
    int *del       = NULL;
    int *add       = NULL;
    int *hit       = NULL;
    int newblossom = 0;
    graph_clique *clique;

    if (depth == 10)
        goto CLEANUP;

    /* Clean up intersecting teeth */

    hit = malloc(graph->nv * sizeof(int));
    del = malloc(graph->nv * sizeof(int));
    add = malloc(graph->nv * sizeof(int));
    check_rval(!hit || !del || !add, "out of memory in work_blossom", CLEANUP);
    newhandle = malloc(graph->nv * sizeof(graph_vertex *));
    check_null(newhandle, "out of memory in work_blossom", CLEANUP);
    newteeth = malloc(graph->na * sizeof(graph_arc *));
    check_null(newteeth, "out of memory in work_blossom", CLEANUP);

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

    if (newblossom)
    {
        rval = grow_ghteeth(graph, repo, newhcount, newhandle, &newtcount,
                            newteeth, yval, cutcount, cuts);
        check_rval(rval, "grow_ghteeth failed", CLEANUP);
    }
    else
    {
        if (hcount >= 3 && (tcount % 2 == 1) && tcount >= 3)
        {
            rval =
            add_blossom(graph, hcount, handle, tcount, teeth, cuts, cutcount);
            check_rval(rval, "add_blossom failed", CLEANUP);
        }
        else if (tcount < 3 && hcount > 1)
        {

            clique = clique_conv_vertices2clique(graph, handle, hcount);
            clique_register_repo(graph, repo, clique);
            clique_free(&clique);
            free(clique);
        }
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
