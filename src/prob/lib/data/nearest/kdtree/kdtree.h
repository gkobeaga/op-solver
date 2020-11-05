#ifndef DATA_KDTREE_H
#define DATA_KDTREE_H

#include "../../data.h"

typedef struct data_kdtree_bnds
{
    double x[2];
    double y[2];
} data_kdtree_bnds;

typedef struct data_kdtree_node data_kdtree_node;
typedef struct data_kdtree_node
{
    double cutval;
    data_kdtree_node *child_lo;
    data_kdtree_node *child_hi;
    data_kdtree_node *parent;
    data_kdtree_bnds *bnds;
    int lo;
    int hi;
    char bucket;
    char empty;
    int cutdim;
} data_kdtree_node;

typedef struct solver_data solver_data;
typedef struct data_kdtree
{
    int n; // Number of nodes
    data_kdtree_node *root;
    data_kdtree_node **nodes;
    int *perm;
    solver_data *data;
} data_kdtree;

void
kdtree_free(data_kdtree **kdtree),
kdtree_delete(data_kdtree *kdtree, int k),
kdtree_delete_all(data_kdtree *kdtree),
kdtree_undelete(data_kdtree *kdtree, int k),
kdtree_undelete_all(data_kdtree *kdtree);

data_kdtree *
kdtree_create(solver_data *data);

#endif
