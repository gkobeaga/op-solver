#include "cp/cp.h"

#define SOLVER_CP_MAX_PATH_CUTS 500

static int
search(solver_graph *graph, graph_vertex *v, int *last, graph_vertex **heap,
       double *len, double *viol, int *npaths, int *pathlen, int **paths);
static int
duplicated(int last, graph_vertex **heap, int npaths, int *pathlen,
           int **paths);
static int
create_path_cut(cp_prob *cp, cp_exact_bac_env *bac_env, cp_cut **cuts,
                int *cutcount, int *verts, int nverts, double *dist2depot);

int
cp_sep_path(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
            cp_cut **cuts)
{
    int rval = 0;
    int i;
    int last;
    graph_vertex **heap;
    double viol, len;
    int **paths;
    int *pathlen;
    int npaths = 0;
    double *dist2depot;

    *cutcount = 0;
    *cuts     = NULL;

    solver_graph *graph = bac_env->ip->lp->sol->graph;

    dist2depot = malloc(graph->nv * sizeof(double));
    memset(dist2depot, 0.0, graph->nv * sizeof(double));
    graph->marker++;
    for (i = 0; i < graph->nv; i++)
    {
        graph->v[i]->used = 0;
        graph->v[i]->mark = graph->marker;
        dist2depot[i]     = data_get_norm(cp->data, 0, i);
    }

    heap    = malloc(graph->nv * sizeof(graph_vertex *));
    paths   = malloc(SOLVER_CP_MAX_PATH_CUTS * sizeof(int *));
    pathlen = malloc(SOLVER_CP_MAX_PATH_CUTS * sizeof(int));
    memset(pathlen, 0, SOLVER_CP_MAX_PATH_CUTS * sizeof(int));

    for (i = 0; i < SOLVER_CP_MAX_PATH_CUTS; i++)
        paths[i] = malloc(graph->nv * sizeof(int));

    for (i = 1; i < graph->nv; i++)
    {
        if (graph->v[i]->y > SOLVER_ZEROPLUS)
        {
            viol = graph->v[i]->y;
            len  = -cp->data->cap;
            last = 0;
            if (search(graph, graph->v[i], &last, heap, &len, &viol, &npaths,
                       pathlen, paths))
                break;
        }
    }

    for (i = 0; i < npaths; i++)
    {
        create_path_cut(cp, bac_env, cuts, cutcount, paths[i], pathlen[i],
                        dist2depot);
    }

    free(dist2depot);
    free(heap);
    for (i = 0; i < SOLVER_CP_MAX_PATH_CUTS; i++) free(paths[i]);
    free(paths);
    free(pathlen);

    return rval;
}

static int
search(solver_graph *graph, graph_vertex *v, int *last, graph_vertex **heap,
       double *len, double *viol, int *npaths, int *pathlen, int **paths)
{
    int i;
    graph_arc *e;
    graph_vertex *other;
    double rhs;

    heap[(*last)++] = v;
    if (*len > 0)
    {
        if (!duplicated(*last - 1, heap, *npaths, pathlen, paths))
        {
            for (i = 0; i < *last - 1; i++) paths[*npaths][i] = heap[i]->i;

            pathlen[*npaths] = *last - 1;
            (*npaths)++;

            if (*npaths == SOLVER_CP_MAX_PATH_CUTS)
            {
                printf(
                " Maximal number %d of required long path was achieved\n",
                SOLVER_CP_MAX_PATH_CUTS);
                return (1);
            }
        }
        (*last)--;
        return (0);
    }

    v->used = 1;
    rhs     = 0.0;
    for (e = v->edge; e; e = outnext(e, v))
    {
        other = otherend(e, v);
        if (!other->used && *len + e->cost <= 0)
            rhs += e->x;
    }

    if (graph->v[0]->used == 0 && *viol - rhs > SOLVER_IP_BAC_MIN_VIOL && *last > 2)
    {
        if (!duplicated(*last, heap, *npaths, pathlen, paths))
        {
            for (i = 0; i < *last; i++) paths[*npaths][i] = heap[i]->i;

            pathlen[*npaths] = *last;
            (*npaths)++;

            if (*npaths == SOLVER_CP_MAX_PATH_CUTS)
            {
                printf(
                " Maximal number %d of required long path was achieved\n",
                SOLVER_CP_MAX_PATH_CUTS);
                return (1);
            }
        }
    }

    *viol -= v->y;

    for (e = v->edge; e; e = outnext(e, v))
    {
        other = otherend(e, v);
        if (other->i != 0 && !other->used)
        {
            *viol += e->x;
            *len += e->cost;

            if (*viol > SOLVER_ZEROPLUS && !other->used)
            {
                if (search(graph, other, last, heap, len, viol, npaths, pathlen,
                           paths))
                    return (1);
            }

            *len -= e->cost;
            *viol -= e->x;
        }
    }

    *viol += v->y;
    v->used = 0;

    (*last)--;
    return 0;
}

static int
duplicated(int last, graph_vertex **heap, int npaths, int *pathlen, int **paths)
{
    int i, k;

    for (k = 0; k < npaths; k++)
    {
        if (last == pathlen[k])
        {
            for (i = 0; i < last; i++)
            {
                if (paths[k][i] != heap[i]->i)
                    break;
            }
            if (i == last)
                return 1;

            for (i = 0; i < last; i++)
            {
                if (paths[k][pathlen[k] - i - 1] != heap[i]->i)
                    break;
            }
            if (i == last)
                return 1;
        }
    }
    return 0;
}

static int
create_path_cut(cp_prob *cp, cp_exact_bac_env *bac_env, cp_cut **cuts,
                int *cutcount, int *verts, int nverts, double *dist2depot)
{
    int rval = 0;
    int i, k;
    double d, cost;
    int depot_in_path = 0;
    graph_arc *arc;
    cp_cut *cut;
    cp_cut_path *path;

    solver_graph *graph = bac_env->ip->lp->sol->graph;

    graph->marker++;

    for (i = 0; i < nverts; i++)
    {
        graph->v[verts[i]]->mark = graph->marker;
        if (verts[i] == 0)
        {
            printf("WARNING: depot inside the path!\n");
            depot_in_path = 1;
        }
    }

    d = -cp->data->cap;
    for (i = 0; i < nverts - 1; i++)
    {
        arc = graph_find_arc_hash(graph->archash, verts[i], verts[i + 1]);
        d += arc->cost;
    }

    k = 0;
    if (depot_in_path)
    {
        if (d <= 0)
        {
            printf("ERROR: depot in a feasible path!\n");
            exit(1);
        }
    }
    else
    {
        for (i = 0; i < graph->nv; i++)
        {
            if (graph->v[i]->mark != graph->marker)
            {
                cost = data_get_norm(cp->data, verts[nverts - 1], i);
                if (dist2depot[verts[0]] + d + cost + dist2depot[i] <= 0)
                    k++;
            }
        }
    }

    cut        = cp_create_cut();
    cut->sense = 'L';

    path        = malloc(sizeof(cp_cut_path));
    path->farcs = NULL;
    path->arcs  = malloc(2 * (nverts - 1) * sizeof(int));

    path->na = nverts - 1;
    for (i = 0; i < nverts - 1; i++)
    {
        path->arcs[2 * i]     = verts[i];
        path->arcs[2 * i + 1] = verts[i + 1];
    }

    path->fna = k;
    if (k)
    {
        path->farcs = malloc(2 * k * sizeof(int));
        for (i = 0, k = 0; i < graph->nv; i++)
        {
            if (graph->v[i]->mark != graph->marker)
            {
                cost = data_get_norm(cp->data, verts[nverts - 1], i);
                if (dist2depot[verts[0]] + d + cost + dist2depot[i] <= 0)
                {
                    path->farcs[2 * k]     = verts[nverts - 1];
                    path->farcs[2 * k + 1] = i;
                    k++;
                }
            }
        }
    }
    cut->path = path;

    cut->next = *cuts;
    *cuts     = cut;
    (*cutcount)++;

    return rval;
}
