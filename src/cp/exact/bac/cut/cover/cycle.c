#include "ip/exact/bac/bac.h"
#include "cp/cp.h"
#include "cp/exact/bac/bac.h"

static void
get_path_from_orig(solver_graph *graph, int target, int comp, double *len,
                   int *na, graph_arc **path);
static int
step_enumerate(int *na, graph_arc **path, double *len, graph_vertex *v,
               graph_arc *e, int target);

int
cp_sep_cover_cycle(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                   cp_cut **cuts)
{
    int rval = 0;
    int i, j;
    solver_graph *tree;
    graph_arc *arc;
    int na1, na2;
    graph_arc **path1 = NULL, **path2 = NULL;
    double len1, len2;

    solver_graph *graph = bac_env->ip->lp->sol->graph;

    cp_cut *cut;
    cp_cut_cover_edge *cover;

    *cutcount = 0;
    *cuts     = NULL;

    tree = graph_get_mst_max(graph);

    path1 = malloc(graph->nv * sizeof(graph_arc *));
    check_null(path1, "out of memory", CLEANUP);
    path2 = malloc(graph->nv * sizeof(graph_arc *));
    check_null(path2, "out of memory", CLEANUP);

    tree->data = cp->data;

    for (i = 0; i < graph->na; i++)
    {
        arc = graph->arcs[i];
        if (!(arc->used) && arc->tail->comp != arc->head->comp)
        {
            if (arc->tail->i)
                get_path_from_orig(tree, arc->tail->i, arc->tail->comp, &len1,
                                   &na1, path1);
            else
            {
                na1  = 0;
                len1 = 0.0;
            }
            if (arc->head->i)
                get_path_from_orig(tree, arc->head->i, arc->head->comp, &len2,
                                   &na2, path2);
            else
            {
                na2  = 0;
                len2 = 0.0;
            }

            if (len1 + len2 + arc->cost > cp->data->cap)
            {

                cut = cp_create_cut();
                check_null(cut, "out of memory", CLEANUP)

                cover       = malloc(sizeof(cp_cut_cover_edge));
                cover->arcs = malloc(2 * (na1 + na2 + 1) * sizeof(int));
                double x    = arc->x;
                for (j = 0; j < na1; j++)
                {
                    cover->arcs[2 * j]     = path1[j]->orig->tail->i;
                    cover->arcs[2 * j + 1] = path1[j]->orig->head->i;
                    x += path1[j]->x;
                }
                for (j = 0; j < na2; j++)
                {
                    cover->arcs[2 * na1 + 2 * j]     = path2[j]->orig->tail->i;
                    cover->arcs[2 * na1 + 2 * j + 1] = path2[j]->orig->head->i;
                    x += path2[j]->x;
                }
                cover->arcs[2 * (na1 + na2)]     = arc->tail->i;
                cover->arcs[2 * (na1 + na2) + 1] = arc->head->i;
                cover->na                        = na1 + na2 + 1;
                cut->sense                       = 'L';
                cut->rhs                         = -1.0;
                cover->strong                    = 1;
                cut->cover_edge                  = cover;

                cut->next = *cuts;
                *cuts     = cut;
                (*cutcount)++;
            }
        }
    }

CLEANUP:
    graph_free(&tree);
    free(path1);
    free(path2);

    return rval;
}

static void
get_path_from_orig(solver_graph *graph, int target, int comp, double *len,
                   int *na, graph_arc **path)
{
    graph_vertex *v;
    graph_arc *e;

    *len = 0.0;
    *na  = 0;

    for (e = graph->v[0]->edge; e; e = outnext(e, graph->v[0]))
    {
        v = otherend(e, graph->v[0]);
        if (v->comp == comp)
            break;
    }

    step_enumerate(na, path, len, graph->v[0], e, target);
}

static int
step_enumerate(int *na, graph_arc **path, double *len, graph_vertex *v,
               graph_arc *e, int target)
{
    int rval;
    graph_arc *f;
    graph_vertex *next, *other;

    path[(*na)++] = e;
    *len += e->cost;

    next = otherend(e, v);
    for (f = next->edge; f; f = outnext(f, next))
    {
        other = otherend(f, next);
        if (other->i != v->i)
        {
            if (other->i == target)
            {
                path[(*na)++] = f;
                *len += f->cost;
                return 1;
            }

            rval = step_enumerate(na, path, len, next, f, target);
            if (rval)
                return 1;
        }
    }

    *len -= e->cost;
    path[--(*na)] = NULL;

    return 0;
}
