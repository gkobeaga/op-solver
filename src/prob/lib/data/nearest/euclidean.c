#include "op-solver.h"

struct shortedge
{
    int end;
    int length;
};

struct table_item
{
    int i;
    int j;
    struct table_item *next;
};

#define dtrunc(x) (((x) > 0.0) ? floor(x) : ceil(x))

static void
node_k_nearest_(solver_data *data, data_kdtree_node *current, int node, int k,
                struct shortedge *nearlist, int *worst_on_list);
static int
ball_in_bounds(solver_data *data, data_kdtree_bnds *bnds, int n, double dist),
put_in_table(struct table_item **table, int i, int j, int *added);

int
data_get_k_nearest_euclidean(solver_data *data, int k)
{
    int rval = 0;
    int i, j, node;
    int added, ntotal = 0;
    int *list                 = NULL;
    struct table_item **table = NULL;
    struct table_item *entry, *next;
    data_map *map = data->map;

    // TODO: skip if already calculated
    map->kn_k      = 0;
    map->kn_ecount = 0;
    map->kn_elist  = NULL;

    table = malloc(data->map->img_n * sizeof(struct table_item *));
    check_null(table, "out of memory", CLEANUP);

    for (i = 0; i < data->map->img_n; i++)
    {
        table[i] = NULL;
    }

    for (node = 0; node < data->map->img_n; node++)
    {
        rval = data_get_node_k_nearest_euclidean(data, node, k, &list);
        check_rval(rval, "data_get_node_k_nearest_euclidean failed", CLEANUP);
        for (i = 0; i < k; i++)
        {
            if (list[i] != -1)
            {
                if (put_in_table(table, node, list[i], &added))
                {
                    fprintf(stderr, "put_in_table failed\n");
                    rval = 1;
                    goto CLEANUP;
                }
                else
                {
                    ntotal += added;
                }
            }
        }
        free(list);
    }

    map->kn_k      = k;
    map->kn_ecount = ntotal;
    map->kn_elist  = malloc(2 * ntotal * sizeof(int));
    check_null(map->kn_elist, "out of memory", CLEANUP);

    for (i = 0, j = 0; i < data->map->img_n; i++)
    {
        for (entry = table[i]; entry; entry = next)
        {
            next                 = entry->next;
            (map->kn_elist)[j++] = i;
            (map->kn_elist)[j++] = entry->j;
            free(entry);
        }
    }

CLEANUP:

    if (table)
        free(table);

    return rval;
    return rval;
}

int
data_get_node_k_nearest_euclidean(solver_data *data, int node, int k,
                                  int **list)
{
    int rval = 0;
    int i;
    data_kdtree_node *current, *last_parent;
    double diff;
    struct shortedge *nearlist = NULL;
    int worst_on_list          = SOLVER_MAXINT;

    nearlist = malloc((k + 1) * sizeof(struct shortedge));
    check_null(nearlist, "", CLEANUP);

    for (i = 0; i < k; i++) nearlist[i].length = SOLVER_MAXINT;
    nearlist[k].length = -SOLVER_MAXINT;

    current = data->map->kdtree->nodes[node];
    node_k_nearest_(data, current, node, k, nearlist, &worst_on_list);
    while (1)
    {
        last_parent = current;
        current     = current->parent;
        if (current == NULL)
            break;
        switch (current->cutdim)
        {
        case 0:
            diff = current->cutval - data->map->x[node];
            if (last_parent == current->child_lo)
            {
                if (worst_on_list > dtrunc(diff))
                    node_k_nearest_(data, current->child_hi, node, k, nearlist,
                                    &worst_on_list);
            }
            else
            {
                if (worst_on_list > dtrunc(-diff))
                    node_k_nearest_(data, current->child_lo, node, k, nearlist,
                                    &worst_on_list);
            }
            break;
        case 1:
            diff = current->cutval - data->map->y[node];
            if (last_parent == current->child_lo)
            {
                if (worst_on_list > dtrunc(diff))
                    node_k_nearest_(data, current->child_hi, node, k, nearlist,
                                    &worst_on_list);
            }
            else
            {
                if (worst_on_list > dtrunc(-diff))
                    node_k_nearest_(data, current->child_lo, node, k, nearlist,
                                    &worst_on_list);
            }
            break;
        }
        if (current->bnds &&
            ball_in_bounds(data, current->bnds, node, worst_on_list))
            break;
    }

    *list = malloc(k * sizeof(int));
    check_null(*list, "out of memory", CLEANUP);

    int ntot = 0;
    for (i = 0; i < k; i++)
    {
        if (nearlist[i].length < SOLVER_MAXINT)
            (*list)[ntot++] = nearlist[i].end;
    }

    if (ntot < k)
    {
        fprintf(stderr, "WARNING: There do not exist %d neighbors\n", k);
        for (i = ntot; i < k; i++) *list[i] = -1;
    }

CLEANUP:
    if (nearlist)
        free(nearlist);
    return rval;
}

static void
node_k_nearest_(solver_data *data, data_kdtree_node *current, int node, int k,
                struct shortedge *nearlist, int *worst_on_list)
{
    int i, j;
    double val, thisx;
    int thisdist;
    data_kdtree *kdtree = data->map->kdtree;

    if (current->bucket)
    {
        for (i = current->lo; i <= current->hi; i++)
        {
            if (kdtree->perm[i] != node)
            {
                thisdist = data_get_norm(data, kdtree->perm[i], node);
                if (*worst_on_list > thisdist)
                // if (*worst_on_list >= thisdist)
                {

                    for (j = 0; nearlist[j + 1].length > thisdist; j++)
                    // for (j = 0; nearlist[j + 1].length >= thisdist; j++)
                    {
                        nearlist[j].end    = nearlist[j + 1].end;
                        nearlist[j].length = nearlist[j + 1].length;
                    }
                    nearlist[j].length = thisdist;
                    nearlist[j].end    = kdtree->perm[i];
                    *worst_on_list     = nearlist[0].length;
                }
            }
        }
    }
    else
    {
        val = current->cutval;
        switch (current->cutdim)
        {
        case 0:
            thisx = data->map->x[node];
            if (thisx < val)
            {
                node_k_nearest_(data, current->child_lo, node, k, nearlist,
                                worst_on_list);
                /* Truncation for floating point coords */
                if (*worst_on_list > dtrunc(val - thisx))
                    node_k_nearest_(data, current->child_hi, node, k, nearlist,
                                    worst_on_list);
            }
            else
            {
                node_k_nearest_(data, current->child_hi, node, k, nearlist,
                                worst_on_list);
                if (*worst_on_list > dtrunc(thisx - val))
                    node_k_nearest_(data, current->child_lo, node, k, nearlist,
                                    worst_on_list);
            }
            break;
        case 1:
            thisx = data->map->y[node];
            if (thisx < val)
            {
                node_k_nearest_(data, current->child_lo, node, k, nearlist,
                                worst_on_list);
                if (*worst_on_list > dtrunc(val - thisx))
                    node_k_nearest_(data, current->child_hi, node, k, nearlist,
                                    worst_on_list);
            }
            else
            {
                node_k_nearest_(data, current->child_hi, node, k, nearlist,
                                worst_on_list);
                if (*worst_on_list > dtrunc(thisx - val))
                    node_k_nearest_(data, current->child_lo, node, k, nearlist,
                                    worst_on_list);
            }
            break;
        }
    }
}

static int
ball_in_bounds(solver_data *data, data_kdtree_bnds *bnds, int n, double dist)
{
    if (dtrunc(data->map->x[n] - bnds->x[0]) < dist ||
        dtrunc(bnds->x[1] - data->map->x[n]) < dist ||
        dtrunc(data->map->y[n] - bnds->y[0]) < dist ||
        dtrunc(bnds->y[1] - data->map->y[n]) < dist)
        return 0;
    return 1;
}
static int
put_in_table(struct table_item **table, int i, int j, int *added)
{
    struct table_item *entry;
    if (j < i)
    {
        int temp;
        SWAP(i, j, temp);
    }

    for (entry = table[i]; entry; entry = entry->next)
    {
        if (entry->j == j)
        {
            *added = 0;
            return 0;
        }
    }

    entry = malloc(sizeof(struct table_item));
    if (!entry)
    {
        *added = 0;
        return 1;
    }

    entry->j    = j;
    entry->next = table[i];
    table[i]    = entry;

    *added = 1;
    return 0;
}
