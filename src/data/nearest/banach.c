#include "op-solver.h"

#define dtrunc(x) (((x) > 0.0) ? floor(x) : ceil(x))

struct shortedge
{
    int length;
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
put_in_table(struct table_item **table, int i, int j, int *added);

int
data_get_k_nearest_banach(solver_data *data, int k)
{
    int rval = 0;
    int i, j, node;
    int added, ntotal = 0;
    int *list                 = NULL;
    struct table_item **table = NULL;
    struct table_item *entry, *next;
    data_map *map = data->map;

    map->kn_k      = 0;
    map->kn_ecount = 0;
    map->kn_elist  = NULL;

    table = malloc(data->n * sizeof(struct table_item *));
    check_null(table, "out of memory", CLEANUP);

    for (i = 0; i < data->n; i++)
    {
        table[i] = NULL;
    }

    for (node = 0; node < data->map->img_n; node++)
    {
        rval = data_get_node_k_nearest_banach(data, node, k, &list);
        check_rval(rval, "failed", CLEANUP);
        for (i = 0; i < k; i++)
        {
            if (list[i] != -1)
            {
                if (put_in_table(table, node, list[i], &added))
                {
                    fprintf(stderr, "failed\n");
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

    if (list)
        free(list);
    if (table)
        free(table);

    return rval;
}

int
data_get_node_k_nearest_banach(solver_data *data, int node, int k, int **list)
{
    int rval = 0;
    int i, j, ntotal;
    struct shortedge *nearlist = NULL;
    int n                      = node;

    if (*list)
    {
        free(*list);
        *list = NULL;
    }

    nearlist = malloc((k + 1) * sizeof(struct shortedge));
    if (!nearlist)
        return 1;
    for (i = 0; i < k; i++) nearlist[i].length = SOLVER_MAXINT;
    nearlist[k].length = -SOLVER_MAXINT;

    for (j = n - 1; j >= 0; --j) list_add(data, n, j, nearlist);

    for (j = n + 1; j < data->map->img_n; j++) list_add(data, n, j, nearlist);

    *list = malloc(k * sizeof(int));
    check_null(*list, "out of memory", CLEANUP);

    ntotal = 0;
    for (i = 0; i < k; i++)
    {
        if (nearlist[i].length < SOLVER_MAXDOUBLE)
            (*list)[ntotal++] = nearlist[i].end;
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
