#include "op-solver.h"

#define CUTOFF 5
#define BNDS_DEPTH 5 /* When bnds info is recorded */

static int
find_max_spread(data_kdtree *kdtree, int lo, int hi);

static data_kdtree_node *
kdtree_create_(data_kdtree *kdtree, int lo, int hi, int depth,
               data_kdtree_bnds *bnds);
static data_kdtree_node *
node_create_(void);
static void
node_free_(data_kdtree_node **node);

data_kdtree *
kdtree_create(solver_data *data)
{
    int rval = 0;
    int i;
    int depth;

    data_kdtree *kdtree    = NULL;
    data_kdtree_bnds *bnds = NULL;

    kdtree = malloc(sizeof(data_kdtree));
    check_null(kdtree, "out of memory", CLEANUP2);
    kdtree->perm  = NULL;
    kdtree->nodes = NULL;
    kdtree->root  = NULL;

    kdtree->n    = data->map->img_n;
    kdtree->data = data;

    kdtree->perm = malloc(kdtree->n * sizeof(int));
    check_null(kdtree->perm, "out of memory", CLEANUP1);
    kdtree->nodes = malloc(kdtree->n * sizeof(data_kdtree_node));
    check_null(kdtree->nodes, "out of memory", CLEANUP1);

    for (i = 0; i < kdtree->n; i++)
    {
        kdtree->perm[i]  = i;
        kdtree->nodes[i] = NULL;
    }

    depth = 0;
    bnds  = malloc(sizeof(data_kdtree_bnds));
    check_null(bnds, "out of memory", CLEANUP1);
    bnds->x[0] = -SOLVER_MAXDOUBLE;
    bnds->x[1] = SOLVER_MAXDOUBLE;
    bnds->y[0] = -SOLVER_MAXDOUBLE;
    bnds->y[1] = SOLVER_MAXDOUBLE;

    kdtree->root = kdtree_create_(kdtree, 0, kdtree->n - 1, depth, bnds);
    check_null(kdtree->root, "kdtree kdtree_create_ failed", CLEANUP1);

    kdtree->root->parent = NULL;

    if (bnds)
        free(bnds);

    return kdtree;

CLEANUP1:
    if (bnds)
        free(bnds);
    if (kdtree->perm)
        free(kdtree->perm);
    if (kdtree->nodes)
        free(kdtree->nodes);
    if (kdtree)
        free(kdtree);
CLEANUP2:
    if (!rval)
        exit(1);
    return NULL;
}

void
kdtree_free(data_kdtree **kdtree)
{
    if (*kdtree)
    {

        node_free_(&((*kdtree)->root));
        if ((*kdtree)->perm)
            free((*kdtree)->perm);

        if ((*kdtree)->nodes)
            free((*kdtree)->nodes);
        free(*kdtree);
        (*kdtree) = NULL;
    }
}

static data_kdtree_node *
node_create_(void)
{
    int rval               = 0;
    data_kdtree_node *node = NULL;

    node = malloc(sizeof(data_kdtree_node));
    check_null(node, "out of memory", CLEANUP);

    node->cutval   = -1;
    node->child_lo = NULL;
    node->child_hi = NULL;
    node->parent   = NULL;
    node->bnds     = NULL;
    node->lo       = -1;
    node->hi       = -1;
    node->bucket   = -1;
    node->empty    = -1;
    node->cutdim   = -1;

    return node;

CLEANUP:
    if (!rval)
        exit(1);
    return NULL;
}

static void
node_free_(data_kdtree_node **node)
{

    if (*node)
    {
        node_free_(&(*node)->child_lo);
        node_free_(&(*node)->child_hi);
        if ((*node)->bnds)
        {
            free((*node)->bnds);
        }
        free(*node);
        *node = NULL;
    }
}

static data_kdtree_bnds *
create_bnd()
{
    int rval               = 0;
    data_kdtree_bnds *bnds = NULL;

    bnds = malloc(sizeof(data_kdtree_bnds));
    check_null(bnds, "out of memory", CLEANUP);

    bnds->x[0] = -1;
    bnds->x[1] = -1;
    bnds->y[0] = -1;
    bnds->y[1] = -1;

    return bnds;

CLEANUP:
    if (!rval)
        exit(1);
    return NULL;
}

static data_kdtree_node *
kdtree_create_(data_kdtree *kdtree, int lo, int hi, int depth,
               data_kdtree_bnds *bnds)
{
    int rval = 0;
    data_kdtree_node *node;
    int i, mid;
    double savebnd;
    data_map *map = kdtree->data->map;

    depth++;

    node = node_create_();
    check_null(node, "kdtree_create_ failed", CLEANUP);

    node->empty = 0;

    if (hi - lo + 1 < CUTOFF)
    {
        node->bucket = 1;
        node->lo     = lo;
        node->hi     = hi;
        for (i = lo; i <= hi; i++) kdtree->nodes[kdtree->perm[i]] = node;
        node->bnds = NULL;
    }
    else
    {
        node->bucket = 0;
        if (!(depth % BNDS_DEPTH))
        {
            node->bnds = create_bnd();
            check_null(node->bnds, "kdtree_create_ failed", CLEANUP);

            node->bnds->x[0] = bnds->x[0];
            node->bnds->x[1] = bnds->x[1];
            node->bnds->y[0] = bnds->y[0];
            node->bnds->y[1] = bnds->y[1];
        }
        else
        {
            node->bnds = NULL;
        }

        node->cutdim = find_max_spread(kdtree, lo, hi);
        mid          = (lo + hi) / 2;
        switch (node->cutdim)
        {
        case 0:
            sort_partial(kdtree->perm, lo, hi, mid, map->x);
            node->cutval = map->x[kdtree->perm[mid]];

            savebnd        = bnds->x[1];
            bnds->x[1]     = node->cutval;
            node->child_lo = kdtree_create_(kdtree, lo, mid, depth, bnds);
            check_null(node->child_lo, "kdtree_create_ failed", CLEANUP);
            bnds->x[1] = savebnd;

            savebnd        = bnds->x[0];
            bnds->x[0]     = node->cutval;
            node->child_hi = kdtree_create_(kdtree, mid + 1, hi, depth, bnds);
            check_null(node->child_hi, "kdtree_create_ failed", CLEANUP);

            bnds->x[0] = savebnd;

            break;

        case 1:
            sort_partial(kdtree->perm, lo, hi, mid, map->y);
            node->cutval = map->y[kdtree->perm[mid]];

            savebnd        = bnds->y[1];
            bnds->y[1]     = node->cutval;
            node->child_lo = kdtree_create_(kdtree, lo, mid, depth, bnds);
            check_null(node->child_lo, "kdtree_create_ failed", CLEANUP);
            bnds->y[1] = savebnd;

            savebnd        = bnds->y[0];
            bnds->y[0]     = node->cutval;
            node->child_hi = kdtree_create_(kdtree, mid + 1, hi, depth, bnds);
            check_null(node->child_hi, "kdtree_create_ failed", CLEANUP);
            bnds->y[0] = savebnd;

            break;
        }
        node->child_lo->parent = node;
        node->child_hi->parent = node;
    }
    return node;
CLEANUP:
    if (!rval)
        exit(1);
    return NULL;
}

static int
find_max_spread(data_kdtree *kdtree, int lo, int hi)
{
    int i;
    double xmax, xmin, xval, xspread;
    double ymax, ymin, yval, yspread;
    data_map *map = kdtree->data->map;

    xmin = map->x[kdtree->perm[lo]];
    xmax = xmin;
    ymin = map->y[kdtree->perm[lo]];
    ymax = ymin;
    for (i = lo + 1; i <= hi; i++)
    {
        xval = map->x[kdtree->perm[i]];
        if (xval < xmin)
            xmin = xval;
        else if (xval > xmax)
            xmax = xval;
        yval = map->y[kdtree->perm[i]];
        if (yval < ymin)
            ymin = yval;
        else if (yval > ymax)
            ymax = yval;
    }

    xspread = xmax - xmin;
    yspread = ymax - ymin;

    if (xspread >= yspread)
        return 0;
    else
        return 1;
}

void
kdtree_delete(data_kdtree *kdtree, int k)
{
    int j, temp;
    data_kdtree_node *node;

    node = kdtree->nodes[k];
    j    = node->lo;
    while (kdtree->perm[j] != k) j++;
    SWAP(kdtree->perm[j], kdtree->perm[node->hi], temp);
    (node->hi)--;
    if (node->lo > node->hi)
    {
        node->empty = 1;
        while ((node = node->parent) && node->child_lo->empty &&
               node->child_hi->empty)
            node->empty = 1;
    }
}

void
kdtree_delete_all(data_kdtree *kdtree)
{
    int i;

    for (i = 0; i < kdtree->n; i++) kdtree_delete(kdtree, i);
}

void
kdtree_undelete(data_kdtree *kdtree, int k)
{
    int j, temp;
    data_kdtree_node *node;

    node = kdtree->nodes[k];
    j    = node->lo;
    while (kdtree->perm[j] != k) j++;
    if (j > node->hi)
    {
        (node->hi)++;
        SWAP(kdtree->perm[j], kdtree->perm[node->hi], temp);
        if (node->empty)
        {
            node->empty = 0;
            while ((node = node->parent) && node->empty) node->empty = 0;
        }
    }
}

void
kdtree_undelete_all(data_kdtree *kdtree)
{
    int i;

    for (i = 0; i < kdtree->n; i++) kdtree_undelete(kdtree, i);
}
