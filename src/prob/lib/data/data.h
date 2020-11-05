#ifndef DATA_H
#define DATA_H

#include "nearest/kdtree/kdtree.h"

/* Mapping */
typedef struct data_map data_map;
typedef struct data_map
{
    int status;
    int img_n;
    int dom_n;
    int *fun;
    int *inv;
    int *orig;
    double *x;
    double *y;
    double *z;
    double *w;
    data_kdtree *kdtree;
    int kn_ecount;
    int *kn_elist;
    int kn_k;
    data_map *prev;
} data_map;

typedef struct solver_data solver_data;
typedef struct solver_graph solver_graph;
typedef struct data_kdtree data_kdtree;
typedef struct solver_data
{
    int n;
    int prob;
    solver_graph *graph;
#define SOLVER_PROB_UNDEFINED (0)
#define SOLVER_PROB_OP (1)
    char name[40];
    double *x;
    double *y;
    double *z;
    double *w;
    int (*edgelen)(solver_data *data, int i, int j);
    data_map *map;
    int from;
    int to;
    double cap;
    int **adj;
    int *adjspace;
    int **len;
    int *lenspace;
    int *degree;
    int *cacheval;
    int *cacheind;
    int cacheM;
    int norm;
    int default_len;   /* for edges not in sparse graph   */
    int sparse_ecount; /* number of edges in sparse graph */
    int nfixed_nodes;
    int nfixed_edges;
    // double   *score;
    double *obj_node;
    double *obj_edge;
    double tot_obj_node;
    double tot_obj_edge;
    int del_ecount;
    int *del_elist;
} solver_data;

#define SOLVER_DATA_TYPE_EUCLIDEAN (128)
#define SOLVER_DATA_TYPE_HILBERT (256)
#define SOLVER_DATA_TYPE_BANACH (512)

#define SOLVER_DATA_TYPE_2D (1024)
#define SOLVER_DATA_TYPE_3D (2048)

#define SOLVER_DATA_NORM_UNDEFINED (0)
#define SOLVER_DATA_NORM_EUCLIDEAN_CEIL                                            \
    (1 | SOLVER_DATA_TYPE_EUCLIDEAN | SOLVER_DATA_TYPE_2D)
#define SOLVER_DATA_NORM_EUCLIDEAN (2 | SOLVER_DATA_TYPE_EUCLIDEAN | SOLVER_DATA_TYPE_2D)
#define SOLVER_DATA_NORM_EUCLIDEAN_3D                                              \
    (3 | SOLVER_DATA_TYPE_HILBERT | /* k-dtree not implemented for 3d */           \
     SOLVER_DATA_TYPE_3D)
#define SOLVER_DATA_NORM_ATT (4 | SOLVER_DATA_TYPE_HILBERT | SOLVER_DATA_TYPE_2D)
#define SOLVER_DATA_NORM_GEOGRAPHIC (5 | SOLVER_DATA_TYPE_HILBERT | SOLVER_DATA_TYPE_2D)
#define SOLVER_DATA_NORM_GEOM (6 | SOLVER_DATA_TYPE_HILBERT | SOLVER_DATA_TYPE_2D)
#define SOLVER_DATA_NORM_MATRIX (7 | SOLVER_DATA_TYPE_BANACH)
#define SOLVER_DATA_NORM_SPARSE (8 | SOLVER_DATA_TYPE_BANACH)

#define SOLVER_DATA_SCALE_GEOGRAPHIC (6378.388 * 3.14 / 180.0)
#define SOLVER_DATA_SCALE_GEOM (6378388.0 * 3.14 / 180.0)
#define SOLVER_DATA_SCALE_ATT (.31622) /* sqrt(1/10) */

solver_data *
data_create(void);
void
data_erase(solver_data *data),
data_free(solver_data **data);

solver_data *
data_read(const char *fname, int format);

int
data_set_norm_type(solver_data *data, int norm),
data_get_norm_type(solver_data *data), data_get_norm(solver_data *data, int i, int j),
data_is_norm_type(solver_data *data, int type);

data_map *
data_create_map(solver_data *data);
data_map *
data_emb_map(solver_data *data, int *vec);
void
data_create_cache(solver_data *data),
data_free_map(data_map **map);

int
data_get_node_k_nearest(solver_data *data, int node, int k, int **list),
data_get_node_k_nearest_euclidean(solver_data *data, int node, int k, int **list),
data_get_node_k_nearest_hilbert(solver_data *data, int node, int k, int **list),
data_get_node_k_nearest_banach(solver_data *data, int node, int k, int **list);

int
data_get_k_nearest(solver_data *data, int k),
data_get_k_nearest_euclidean(solver_data *data, int k),
data_get_k_nearest_hilbert(solver_data *data, int k),
data_get_k_nearest_banach(solver_data *data, int k);
void
data_plot_nearest(solver_data *data);

#endif
