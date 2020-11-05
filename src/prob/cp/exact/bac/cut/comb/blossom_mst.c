#include "cp/cp.h"

#define inhandle(arc, marker)                                                  \
    ((arc)->tail->mark == marker ? (arc)->tail : (arc)->head)
#define outhandle(arc, marker)                                                 \
    ((arc)->tail->mark != marker ? (arc)->tail : (arc)->head)

static int
grab_component(solver_graph *tree, int *hcount, graph_vertex **handle,
               graph_vertex *v),
grow_teeth(solver_graph *graph, int hcount, graph_vertex **handle, int *tcount,
           graph_arc **teeth),
work_blossom(solver_graph *graph, int hcount, graph_vertex **handle, int tcount,
             graph_arc **teeth, int depth, cp_cut **cuts, int *cutcount),
add_blossom(solver_graph *graph, int hcount, graph_vertex **handle, int tcount,
            graph_arc **teeth, cp_cut **cuts, int *cutcount);

static int
sort_arcs(const void *xx, const void *yy)
{
    graph_arc *x = *(graph_arc **)xx, *y = *(graph_arc **)yy;

    if (x->x < y->x)
        return -1;
    if (x->x > y->x)
        return +1;

    return 0;
}

int
cp_sep_blossom_mst(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                   cp_cut **cuts)
{
    int rval = 0;
    int i, j;
    int hcount, tcount;
    graph_vertex **handle   = NULL;
    graph_arc **teeth       = NULL;
    graph_arc **sorted_arcs = NULL;
    graph_clique *clique    = NULL;
    graph_clique_repo *repo = NULL;
    cp_cut_list *cut_list   = NULL;
    graph_arc *e;

    solver_graph *graph = bac_env->ip->lp->sol->graph;

    *cutcount = 0;
    *cuts     = NULL;

    repo = clique_create_repo(graph->nv);
    check_null(repo, "out of memory ", CLEANUP);
    cut_list = cp_create_cut_list(bac_env->param->sec_max_cut, 2.0);
    check_null(cut_list, "failed", CLEANUP);

    for (i = 0; i < graph->nv; i++) graph->v[i]->comp = i;

    sorted_arcs = malloc(graph->na * sizeof(graph_arc *));
    for (i = 0; i < graph->na; i++) sorted_arcs[i] = graph->arcs[i];
    qsort(sorted_arcs, graph->na, sizeof(graph_arc *), sort_arcs);

    handle = malloc(graph->nv * sizeof(graph_vertex *));
    check_null(handle, "out of memory", CLEANUP);
    teeth = malloc(graph->na * sizeof(graph_arc *));
    check_null(teeth, "out of memory", CLEANUP);

    for (i = 0; i < graph->na - 1; i++)
    {
        e = sorted_arcs[i];

        // As suggested in (Fischetti, 1998).
        if (e->x >= SOLVER_ONEMINUS)
            goto OUT;

        if (e->tail->comp != e->head->comp)
        {
            hcount = 0;
            for (j = 0; j < graph->nv; j++)
            {
                if (graph->v[j]->comp == e->head->comp)
                {
                    graph->v[j]->comp = e->tail->comp;
                    handle[hcount++]  = graph->v[j];
                }
                else if (graph->v[j]->comp == e->tail->comp)
                    handle[hcount++] = graph->v[j];
            }

            grow_teeth(graph, hcount, handle, &tcount, teeth);

            if (tcount > 2)
            {
                rval = work_blossom(graph, hcount, handle, tcount, teeth, 0,
                                    cuts, cutcount);
                check_rval(rval, "work_blossom failed", CLEANUP);
            }
            else if (tcount < 3 && hcount > 2)
            {
                clique = clique_conv_vertices2clique(graph, handle, hcount);
                clique_register_repo(graph, repo, clique);
                clique_free(&clique);
                free(clique);
            }
        }
    }
OUT:

    cp_get_clique_repo_sec_cuts(cp, bac_env, graph, repo, cut_list);
    cp_get_cut_list(cp, cut_list, cutcount, cuts);

CLEANUP:
    if (handle)
        free(handle);
    if (teeth)
        free(teeth);
    if (repo)
        clique_free_repo(&repo);
    if (cut_list)
        cp_free_cut_list(&cut_list);
    if (sort_arcs)
        free(sorted_arcs);

    return rval;
}

static int
grow_teeth(solver_graph *graph, int hcount, graph_vertex **handle, int *tcount,
           graph_arc **teeth)
{
    int i;
    graph_vertex *v;
    graph_arc *arc;
    *tcount      = 0;
    int tmpcount = 0;

    graph->marker++;
    for (i = 0; i < hcount; i++) handle[i]->mark = graph->marker;

    if (hcount < 3)
        goto CLEANUP;

    for (i = 0; i < hcount; i++)
    {
        v = handle[i];
        for (arc = v->edge; arc; arc = outnext(arc, v))
        {
            if (arc->tail->mark != arc->head->mark)
                teeth[(tmpcount)++] = arc;
        }
    }

    *tcount = tmpcount;

CLEANUP:

    return 0;
}

static int
work_blossom(solver_graph *graph, int hcount, graph_vertex **handle, int tcount,
             graph_arc **teeth, int depth, cp_cut **cuts, int *cutcount)
{
    graph_vertex *in, *out;
    graph_arc **newteeth     = NULL;
    graph_vertex **newhandle = NULL;
    int i, newhcount, newtcount, rval = 0;
    int *del       = NULL;
    int *add       = NULL;
    int *hit       = NULL;
    int newblossom = 0;

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

    grow_teeth(graph, newhcount, newhandle, &newtcount, newteeth);

    {
        if (newblossom)
        {
            rval = work_blossom(graph, newhcount, newhandle, newtcount,
                                newteeth, depth + 1, cuts, cutcount);
            check_rval(rval, "work_blossom failed", CLEANUP);
        }
        else if (newhcount >= 3 && (newtcount % 2 == 1) && newtcount >= 3)
        {
            rval = add_blossom(graph, newhcount, newhandle, newtcount, newteeth,
                               cuts, cutcount);
            check_rval(rval, "add_blossom failed", CLEANUP);
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
sort_teeth_ind(const void *xx, const void *yy)
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
    check_null(verts, "out of memory ", CLEANUP);

    cut = cp_create_cut();
    check_null(cut, "out of memory ", CLEANUP);

    cut->handles = malloc(sizeof(graph_clique *));
    check_null(cut->handles, "out of memory", CLEANUP);
    cut->teeth = malloc(tcount * sizeof(graph_clique *));

    check_null(cut->teeth, "out of memory", CLEANUP);
    cut->tcount = tcount;
    cut->hcount = 1;

    cut->handles[0] = clique_conv_vertices2clique(graph, handle, hcount);
    qsort(teeth, tcount, sizeof(graph_arc *), sort_teeth_ind);

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
    check_rval(rval, "failed", CLEANUP);

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
