#include "op-solver.h"

int
data_get_k_nearest(solver_data *data, int k)
{
    int rval      = 0;
    data_map *map = data->map;

    if (k >= map->img_n)
    {
        map->kn_ecount = map->img_n * (map->img_n - 1) / 2;
        map->kn_elist  = malloc(2 * (map->kn_ecount) * sizeof(int));
        check_null(map->kn_elist, " failed", CLEANUP);
        for (int i = 0, l = 0; i < map->img_n - 1; i++)
        {
            for (int j = i + 1; j < map->img_n; j++)
            {
                map->kn_elist[2 * l]     = i;
                map->kn_elist[2 * l + 1] = j;
                l++;
            }
        }
        map->kn_k = map->img_n;
        return rval;
    }

    if (data_is_norm_type(data, SOLVER_DATA_TYPE_EUCLIDEAN))
    {
        rval = data_get_k_nearest_euclidean(data, k);
        check_rval(rval, "kdtree_k_nearest failed", CLEANUP);
    }
    else if (data_is_norm_type(data, SOLVER_DATA_TYPE_HILBERT))
    {
        rval = data_get_k_nearest_hilbert(data, k);
        check_rval(rval, "kdtree_k_nearest failed", CLEANUP);
    }
    else if (1 || data_is_norm_type(data, SOLVER_DATA_TYPE_BANACH))
    {
        rval = data_get_k_nearest_banach(data, k);
        check_rval(rval, "kdtree_k_nearest failed", CLEANUP);
    }
    else
    {
        rval = 1;
        fprintf(stderr, "data_get_k_nearest failed\n");
    }

CLEANUP:
    return rval;
}

int
data_get_node_k_nearest(solver_data *data, int node, int k, int **list)
{
    int rval      = 0;
    data_map *map = data->map;

    if (k >= map->img_n)
    {
        *list = malloc(map->img_n * sizeof(int));
        check_null(*list, "kdtree_node_k_nearest failed", CLEANUP);
        for (int i = 0; i < map->img_n; i++)
        {
            if (i != node)
                *list[i] = i;
        }
        return rval;
    }

    if (data_is_norm_type(data, SOLVER_DATA_TYPE_EUCLIDEAN))
    {
        rval = data_get_node_k_nearest_euclidean(data, node, k, list);
        check_rval(rval, "kdtree_k_nearest failed", CLEANUP);
    }
    else if (data_is_norm_type(data, SOLVER_DATA_TYPE_HILBERT))
    {
        rval = data_get_node_k_nearest_hilbert(data, node, k, list);
        check_rval(rval, "kdtree_k_nearest failed", CLEANUP);
    }
    else if (data_is_norm_type(data, SOLVER_DATA_TYPE_BANACH))
    {
        rval = data_get_node_k_nearest_banach(data, node, k, list);
        check_rval(rval, "kdtree_k_nearest failed", CLEANUP);
    }
    else
    {
        rval = 1;
        fprintf(stderr, "data_get_node_k_nearest failed\n");
    }

CLEANUP:
    return rval;
}

void
data_plot_nearest(solver_data *data)
{
    solver_graph *graph = graph_create();
    graph_arc *arc;
    graph_add_vertices(graph, data->n);
    graph->data = data;

    for (int i = 0; i < data->map->kn_ecount; i++)
    {
        arc          = graph_add_arc(graph, data->map->kn_elist[2 * i],
                                     data->map->kn_elist[2 * i + 1]);
        arc->x       = 1;
        arc->tail->y = 1;
        arc->head->y = 1;
    }

    graph_plot(graph);

    graph_free(&graph);
}
