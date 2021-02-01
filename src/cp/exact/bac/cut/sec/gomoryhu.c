#include "cp/cp.h"
#include "cp/exact/bac/bac.h"

static int
test_cuts(cp_prob *cp, solver_graph *graph, graph_ghtree_node *parent,
          graph_clique_repo *repo);
static void
get_verts(graph_ghtree_node *parent, int *count, graph_vertex **verts);

int
cp_sep_sec_exact_gomoryhu(cp_prob *cp, cp_exact_bac_env *bac_env,
                          solver_graph *graph, graph_clique_repo *repo)
{
    int rval = 0;
    graph_ghtree_node *child;

    graph_ghtree *ghtree;

    graph_get_ghtree(graph, graph->tail->i, &ghtree);

    for (child = ghtree->root; child; child = child->next_sibling)
        test_cuts(cp, graph, child, repo);

    graph_free_ghtree(ghtree);

    return rval;
}

static int
test_cuts(cp_prob *cp, solver_graph *graph, graph_ghtree_node *parent,
          graph_clique_repo *repo)
{
    int rval = 0;
    int i;
    graph_clique *clique  = NULL;
    graph_vertex **cverts = NULL;
    graph_vertex **verts  = NULL;
    graph_vertex *v;
    graph_vertex **tverts = NULL;
    int nverts;
    int cvcount = 0;
    graph_ghtree_node *child;
    int fixed_in  = 0;
    double obj_in = 0;

    tverts = malloc(graph->n3v * sizeof(graph_vertex *));
    check_null(tverts, "out of memory", CLEANUP);

    if (parent->cutval - 2 * parent->in_max->y < 0 - SOLVER_ZEROPLUS &&
        parent->ndesc > 1)
    {
        verts = malloc(graph->nv * sizeof(graph_vertex *));
        check_null(verts, "out of memory", CLEANUP);
        cverts = malloc(parent->ndesc * sizeof(graph_vertex *));
        check_null(cverts, "out of memory", CLEANUP);

        cvcount = 0;
        get_verts(parent, &cvcount, cverts);
        assert(cvcount == parent->ndesc);

        for (i = 0, nverts = 0; i < cvcount; i++)
        {
            verts[nverts++] = cverts[i];
            obj_in += cverts[i]->obj;
            if (cverts[i]->fixed)
                fixed_in++;
            for (v = cverts[i]->shrunk->members; v; v = v->members)
            {
                verts[nverts++] = v->orig;
                obj_in += v->obj;
                if (v->fixed)
                    fixed_in++;
            }
        }

        if (nverts <= graph->nv / 2 && nverts > 2 && nverts < graph->nv - 2)
            clique = clique_conv_vertices2clique(graph, verts, nverts);
        else
            clique = clique_conv_vertices2coclique(graph, verts, nverts);
        clique->val = parent->cutval;

        clique_register_repo(graph, repo, clique);

        clique_free(&clique);
        free(verts);
        free(cverts);
    }

    else if (1 && cp->data->tot_obj_edge == 0.0)
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

        cp_find_sec_cliques(graph, tverts, tcount, parent->cutval,
                            cp->ip->lowerboundG, repo);
    }

    for (child = parent->child; child; child = child->next_sibling)
    {
        test_cuts(cp, graph, child, repo);
    }

CLEANUP:
    if (tverts)
        free(tverts);
    return rval;
}

static void
get_verts(graph_ghtree_node *parent, int *count, graph_vertex **verts)
{
    graph_ghtree_node *child;

    verts[(*count)++] = parent->special->orig;

    for (child = parent->child; child; child = child->next_sibling)
    {
        get_verts(child, count, verts);
    }

    return;
}
