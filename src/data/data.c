#include "op-solver.h"

static double
dtrunc(double);

static int
edgelen_nonorm(solver_data *data, int i, int j),
euclid_edgelen(solver_data *data, int i, int j),
euclid_ceiling_edgelen(solver_data *data, int i, int j),
euclid3d_edgelen(solver_data *data, int i, int j),
geographic_edgelen(solver_data *data, int i, int j),
geom_edgelen(solver_data *data, int i, int j),
att_edgelen(solver_data *data, int i, int j),
matrix_edgelen(solver_data *data, int i, int j),
sparse_edgelen(solver_data *data, int i, int j);

static int
edgelen_nonorm(solver_data *data, int i, int j)
{
    if (i != 0 || j != 0 || data)
    {
        fprintf(stderr, "This is a FATAL ERROR\n");
        exit(1);
    }
    return -1;
}

static void
data_init(solver_data *data)
{
    data->n             = 0;
    data->graph         = NULL;
    data->x             = NULL;
    data->y             = NULL;
    data->z             = NULL;
    data->w             = NULL;
    data->prob          = SOLVER_PROB_UNDEFINED;
    data->map           = NULL;
    data->adj           = NULL;
    data->adjspace      = NULL;
    data->len           = NULL;
    data->lenspace      = NULL;
    data->degree        = NULL;
    data->norm          = SOLVER_DATA_NORM_UNDEFINED;
    data->default_len   = 100000;
    data->sparse_ecount = 0;
    data->obj_node      = NULL;
    data->obj_edge      = NULL;
    data->tot_obj_node  = 0.0;
    data->tot_obj_edge  = 0.0;
    data->edgelen       = edgelen_nonorm;
    data->del_elist     = NULL;
    data->del_ecount    = 0;
    data->cacheind      = NULL;
    data->cacheval      = NULL;
    return;
}

solver_data *
data_create(void)
{
    solver_data *data = malloc(sizeof(solver_data));
    data_init(data);
    return data;
}

static void
data_delete(solver_data *data);

void
data_erase_data(solver_data *data)
{
    data_delete(data);
    data_init(data);
    return;
}

static void
data_delete(solver_data *data)
{
    if (data->graph)
        graph_free(&(data->graph));
    if (data->x)
        free(data->x);
    if (data->y)
        free(data->y);
    if (data->z)
        free(data->z);
    if (data->w)
        free(data->w);
    if (data->adj)
        free(data->adj);
    if (data->adjspace)
        free(data->adjspace);
    if (data->len)
        free(data->len);
    if (data->lenspace)
        free(data->lenspace);
    if (data->degree)
        free(data->degree);
    if (data->obj_node)
        free(data->obj_node);
    if (data->obj_edge)
        free(data->obj_edge);

    data->del_ecount = 0;
    if (data->del_elist)
        free(data->del_elist);

    if (data->map)
        data_free_map(&data->map);
    if (data->cacheval)
        free(data->cacheval);
    if (data->cacheind)
        free(data->cacheind);
    return;
}

void
data_free(solver_data **data)
{
    data_delete(*data);
    free(*data);
    *data = NULL;
    return;
}

void
data_create_cache(solver_data *data)
{
    int i;

    if (data->n == 0)
    {
        fprintf(stderr, "null size problem in data_init_cache\n");
        exit(1);
    }

    if (data == NULL)
    {
        fprintf(stderr, "data not initialized in data_init_cache\n");
        exit(1);
    }

    i = 0;
    while ((1 << i) < (data->n << 2)) i++;
    data->cacheM = (1 << i);

    data->cacheind = malloc(data->cacheM * sizeof(int));
    data->cacheval = malloc(data->cacheM * sizeof(int));
    memset(data->cacheind, -1, data->cacheM * sizeof(int));

    if (!data->cacheind || !data->cacheval)
    {
        fprintf(stderr, "out of memory in init_cache\n");
        exit(1);
    }

    data->cacheM--;

    return;
}

data_map *
data_create_map(solver_data *data)
{
    data_map *map;
    data_map *oldmap = data->map;

    map       = malloc(sizeof(data_map));
    map->fun  = malloc(oldmap->img_n * sizeof(int));
    map->inv  = malloc(oldmap->img_n * sizeof(int));
    map->orig = malloc(oldmap->img_n * sizeof(int));
    if (data->x)
        map->x = malloc(oldmap->img_n * sizeof(double));
    else
        map->x = NULL;
    if (data->y)
        map->y = malloc(oldmap->img_n * sizeof(double));
    else
        map->y = NULL;
    if (data->z)
        map->z = malloc(oldmap->img_n * sizeof(double));
    else
        map->z = NULL;
    if (data->w)
        map->w = malloc(oldmap->img_n * sizeof(double));
    else
        map->w = NULL;
    for (int i = 0; i < oldmap->img_n; i++)
    {
        map->fun[i] = map->inv[i] = map->orig[i] = i;
        if (data->x)
            map->x[i] = oldmap->x[i];
        if (data->y)
            map->y[i] = oldmap->y[i];
        if (data->z)
            map->z[i] = oldmap->z[i];
        if (data->w)
            map->w[i] = oldmap->w[i];
    }
    map->dom_n     = oldmap->img_n;
    map->img_n     = oldmap->img_n;
    map->status    = oldmap->status + 1;
    map->prev      = oldmap;
    data->map      = map;
    map->kn_elist  = NULL;
    map->kn_ecount = 0;
    map->kn_k      = 0;
    map->kdtree    = NULL;

    return map;
}

void
data_free_map(data_map **map)
{
    if (*map)
    {
        if ((*map)->fun)
            free((*map)->fun);
        if ((*map)->inv)
            free((*map)->inv);
        if ((*map)->orig)
            free((*map)->orig);
        if ((*map)->x)
            free((*map)->x);
        if ((*map)->y)
            free((*map)->y);
        if ((*map)->z)
            free((*map)->z);
        if ((*map)->w)
            free((*map)->w);
        kdtree_free(&(*map)->kdtree);

        (*map)->kn_ecount = 0;
        if ((*map)->kn_elist)
            free((*map)->kn_elist);
        free(*map);
        *map = NULL;
    }
    return;
}

data_map *
data_emb_map(solver_data *data, int *vec)
{

    data_map *oldmap = data->map;
    data_map *map    = data_create_map(data);
    if (!map)
        return NULL;

    int n = 0;
    for (int i = 0; i < oldmap->img_n; i++)
    {
        if (vec[i])
        {
            map->fun[i]  = n;
            map->inv[n]  = i;
            map->orig[n] = oldmap->orig[i];
            if (data->x)
                map->x[n] = oldmap->x[i];
            if (data->y)
                map->y[n] = oldmap->y[i];
            if (data->z)
                map->z[n] = oldmap->z[i];
            if (data->w)
                map->w[n] = oldmap->w[i];
            n++;
        }
    }

    map->kn_elist  = NULL;
    map->kn_ecount = 0;
    map->kn_k      = 0;
    map->img_n     = n;
    map->dom_n     = oldmap->img_n;
    data->map      = map;
    if (oldmap->kdtree)
        map->kdtree = kdtree_create(data);
    return map;
}

int
data_get_norm(solver_data *data, int i, int j)
{
    int ind;
    int orig_i = data->map->orig[i];
    int orig_j = data->map->orig[j];

    if (orig_i > orig_j)
    {
        int temp;
        SWAP(orig_i, orig_j, temp);
    }

    ind = (((orig_i << 8) + orig_i + orig_j) & (data->cacheM));

    if (data->cacheind[ind] != orig_i)
    {
        data->cacheind[ind] = orig_i;
        data->cacheval[ind] = (data->edgelen)(data, orig_i, orig_j);
    }
    return data->cacheval[ind];
}

int
data_set_norm_type(solver_data *data, int norm)
{
    switch (norm)
    {
    case SOLVER_DATA_NORM_EUCLIDEAN_CEIL:
        data->edgelen = euclid_ceiling_edgelen;
        break;
    case SOLVER_DATA_NORM_EUCLIDEAN:
        data->edgelen = euclid_edgelen;
        break;
    case SOLVER_DATA_NORM_EUCLIDEAN_3D:
        data->edgelen = euclid3d_edgelen;
        break;
    case SOLVER_DATA_NORM_GEOGRAPHIC:
        data->edgelen = geographic_edgelen;
        break;
    case SOLVER_DATA_NORM_GEOM:
        data->edgelen = geom_edgelen;
        break;
    case SOLVER_DATA_NORM_ATT:
        data->edgelen = att_edgelen;
        break;
    case SOLVER_DATA_NORM_MATRIX:
        data->edgelen = matrix_edgelen;
        break;
    case SOLVER_DATA_NORM_SPARSE:
        data->edgelen = sparse_edgelen;
        break;
    default:
        printf("ERROR:  Unknown NORM %d.\n", norm);
        return 1;
    }
    data->norm = norm;

    return 0;
}

int
data_get_norm_type(solver_data *data)
{
    return data->norm;
}

int
data_is_norm_type(solver_data *data, int type)
{
    return (data->norm & type) == type;
}

#if 0
static int
max_edgelen(solver_data *data, int i, int j)
{
    double t1 = data->x[i] - data->x[j], t2 = data->y[i] - data->y[j];

    if (t1 < 0)
        t1 *= -1;
    if (t2 < 0)
        t2 *= -1;
    t1 += 0.5;
    t2 += 0.5;

    return (int)(t1 < t2 ? t2 : t1);
}
#endif

static int
euclid_edgelen(solver_data *data, int i, int j)
{
    double t1 = data->x[i] - data->x[j], t2 = data->y[i] - data->y[j];
    int temp;

    temp = (int)(sqrt(t1 * t1 + t2 * t2) + 0.5);
    return temp;

    return (int)((temp) + 0.5);
}

static int
euclid3d_edgelen(solver_data *data, int i, int j)
{
    double t1 = data->x[i] - data->x[j], t2 = data->y[i] - data->y[j];
    double t3 = data->z[i] - data->z[j];
    int temp;

    temp = (int)(sqrt(t1 * t1 + t2 * t2 + t3 * t3) + 0.5);
    return temp;
}

static int
euclid_ceiling_edgelen(solver_data *data, int i, int j)
{
    double t1 = data->x[i] - data->x[j], t2 = data->y[i] - data->y[j];
    return (int)(ceil(sqrt(t1 * t1 + t2 * t2)));
}

#define GH_PI (3.141592)

static int
geographic_edgelen(solver_data *data, int i, int j)
{
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    int dd;
    double x1 = data->x[i], x2 = data->x[j], yy1 = data->y[i], yy2 = data->y[j];

    deg  = dtrunc(x1);
    min  = x1 - deg;
    lati = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg  = dtrunc(x2);
    min  = x2 - deg;
    latj = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg   = dtrunc(yy1);
    min   = yy1 - deg;
    longi = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg   = dtrunc(yy2);
    min   = yy2 - deg;
    longj = GH_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos(longi - longj);
    q2 = cos(lati - latj);
    q3 = cos(lati + latj);
    dd =
    (int)(6378.388 * acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    return dd;
}

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int
geom_edgelen(solver_data *data, int i, int j)
{
    double lati, latj, longi, longj;
    double q1, q2, q3, q4, q5;

    lati = M_PI * data->x[i] / 180.0;
    latj = M_PI * data->x[j] / 180.0;

    longi = M_PI * data->y[i] / 180.0;
    longj = M_PI * data->y[j] / 180.0;

    q1 = cos(latj) * sin(longi - longj);
    q3 = sin((longi - longj) / 2.0);
    q4 = cos((longi - longj) / 2.0);
    q2 = sin(lati + latj) * q3 * q3 - sin(lati - latj) * q4 * q4;
    q5 = cos(lati - latj) * q4 * q4 - cos(lati + latj) * q3 * q3;
    return (int)(6378388.0 * atan2(sqrt(q1 * q1 + q2 * q2), q5) + 1.0);
}

static int
att_edgelen(solver_data *data, int i, int j)
{
    double xd  = data->x[i] - data->x[j];
    double yd  = data->y[i] - data->y[j];
    double rij = sqrt((xd * xd + yd * yd) / 10.0);
    double tij = dtrunc(rij);
    int dij;

    if (tij < rij)
        dij = (int)tij + 1;
    else
        dij = (int)tij;
    return dij;
}

static double
dtrunc(double x)
{
    int k;

    k = (int)x;
    x = (double)k;
    return x;
}

static int
matrix_edgelen(solver_data *data, int i, int j)
{
    if (i > j)
        return (data->adj[i])[j];
    else
        return (data->adj[j])[i];
}

static int
sparse_edgelen(solver_data *data, int i, int j)
{
    int *adj;
    int k, deg;

    if (i > j)
        SWAP(i, j, k);

    adj = data->adj[i];
    deg = data->degree[i];

    for (k = 0; k < deg; k++)
    {
        if (adj[k] == j)
            return data->len[i][k];
    }
    return data->default_len;
}
