#include "op-solver.h"

#define dtrunc(x) (((x) > 0.0) ? floor(x) : ceil(x))

struct x_sorted_item
{
    int i;
    double x;
};

struct shortedge
{
    double length;
    int end;
};

struct table_item
{
    int i;
    int j;
    struct table_item *next;
};

static void
list_add(solver_data *data, int n, int m, struct shortedge *nearlist);

static int
put_in_table(struct table_item **table, int i, int j, int *added),
sort_x(const void *xx, const void *yy);

int
data_get_k_nearest_hilbert(solver_data *data, int k)
{
    int rval = 0;
    int i, j, node;
    int added, ntotal = 0;
    int *list                      = NULL;
    struct x_sorted_item *x_sorted = NULL;
    struct table_item **table      = NULL;
    data_map *oldmap               = data->map;
    data_map *newmap               = NULL;
    struct table_item *entry, *next;

    // TODO: skip if already calculated
    oldmap->kn_k      = 0;
    oldmap->kn_ecount = 0;
    oldmap->kn_elist  = NULL;

    x_sorted = malloc(data->n * sizeof(struct x_sorted_item));
    check_null(x_sorted, "map_create_hilbert failed", CLEANUP);
    table = malloc(data->n * sizeof(struct table_item *));
    check_null(table, "out of memory", CLEANUP);

    for (i = 0; i < oldmap->img_n; i++)
    {
        x_sorted[i].i = i;
        x_sorted[i].x = oldmap->x[i];
        table[i]      = NULL;
    }

    newmap = data_create_map(data);
    check_null(newmap, "out of memory", CLEANUP);

    qsort(x_sorted, oldmap->img_n, sizeof(struct x_sorted_item), sort_x);

    for (i = 0; i < oldmap->img_n; i++)
    {
        newmap->fun[x_sorted[i].i] = i;
        newmap->inv[i]             = x_sorted[i].i;
        newmap->orig[i]            = oldmap->orig[x_sorted[i].i];
    }

    for (node = 0; node < newmap->dom_n; node++)
    {
        rval = data_get_node_k_nearest_hilbert(data, node, k, &list);
        check_rval(rval, "failed", CLEANUP);
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
    }

    oldmap->kn_k      = k;
    oldmap->kn_ecount = ntotal;
    oldmap->kn_elist  = malloc(2 * ntotal * sizeof(int));
    check_null(oldmap->kn_elist, "out of memory", CLEANUP);

    for (i = 0, j = 0; i < data->map->dom_n; i++)
    {
        for (entry = table[i]; entry; entry = next)
        {
            next                    = entry->next;
            (oldmap->kn_elist)[j++] = i;
            (oldmap->kn_elist)[j++] = entry->j;
            free(entry);
        }
    }

CLEANUP:

    data->map = oldmap;

    if (x_sorted)
        free(x_sorted);
    if (list)
        free(list);
    if (table)
        free(table);
    if (newmap)
        data_free_map(&newmap);

    return rval;
}

int
data_get_node_k_nearest_hilbert(solver_data *data, int node, int k, int **list)
{
    int rval = 0;
    int i, j, ntotal;
    struct shortedge *nearlist = NULL;
    double scale;
    int n = data->map->fun[node];
    int norm;

    if (*list)
    {
        free(*list);
        *list = NULL;
    }

    nearlist = malloc((k + 1) * sizeof(struct shortedge));
    if (!nearlist)
        return 1;
    for (i = 0; i < k; i++) nearlist[i].length = SOLVER_MAXDOUBLE;
    nearlist[k].length = -SOLVER_MAXDOUBLE;

    norm = data_get_norm_type(data);
    if (norm == SOLVER_DATA_NORM_GEOGRAPHIC)
        scale = SOLVER_DATA_SCALE_GEOGRAPHIC;
    else if (norm == SOLVER_DATA_NORM_GEOM)
        scale = SOLVER_DATA_SCALE_GEOM;
    else if (norm == SOLVER_DATA_NORM_ATT)
        scale = SOLVER_DATA_SCALE_ATT;
    else
        scale = 1.0;

    for (j = n - 1; j >= 0 && dtrunc((data->map->x[n] - data->map->x[j]) *
                                     scale) < nearlist[0].length;
         --j)
    {
        list_add(data, n, j, nearlist);
    }

    for (j = n + 1; j < data->map->dom_n &&
                    dtrunc((data->map->x[j] - data->map->x[n]) * scale) <
                    nearlist[0].length;
         j++)
    {
        list_add(data, n, j, nearlist);
    }

    *list = malloc(k * sizeof(int));
    check_null(*list, "out of memory", CLEANUP);

    ntotal = 0;
    for (i = 0; i < k; i++)
    {
        if (nearlist[i].length < SOLVER_MAXDOUBLE)
            (*list)[ntotal++] = data->map->inv[nearlist[i].end];
    }

    if (ntotal < k)
    {
        fprintf(stderr, "WARNING: There do not exist %d neighbors\n", k);
        return 1;
    }

CLEANUP:
    if (nearlist)
        free(nearlist);
    return rval;
}

static void
list_add(solver_data *data, int n, int m, struct shortedge *nearlist)
{
    int i;
    int thisdist;

    thisdist = data_get_norm(data, n, m);

    if (thisdist < nearlist[0].length)
    {
        for (i = 0; nearlist[i + 1].length > thisdist; i++)
        {
            nearlist[i].end    = nearlist[i + 1].end;
            nearlist[i].length = nearlist[i + 1].length;
        }
        nearlist[i].length = thisdist;
        nearlist[i].end    = m;
    }
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

static int
sort_x(const void *xx, const void *yy)
{
    struct x_sorted_item *x = (struct x_sorted_item *)xx,
                         *y = (struct x_sorted_item *)yy;

    if (x->x < y->x)
        return +1;
    if (x->x > y->x)
        return -1;

    return 0;
}
