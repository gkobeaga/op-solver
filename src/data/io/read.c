#include "op-solver.h"
#include "../nearest/kdtree/kdtree.h"

#define MATRIX_LOWER_DIAG_ROW 0
#define MATRIX_UPPER_ROW 1
#define MATRIX_UPPER_DIAG_ROW 2
#define MATRIX_FULL_MATRIX 3

static data_map *
create_first_map(solver_data *data);

static int
data_read_(solver_data *data, FILE *file)
{
    int rval       = 0;
    int matrixform = MATRIX_LOWER_DIAG_ROW;
    char buf[256], key[256], field[256];
    char *p;
    printf("\n");
    while (fgets(buf, 254, file))
    {
        p = buf;
        while (*p != '\0')
        {
            if (*p == ':')
                *p = ' ';
            p++;
        }
        p = buf;
        if (sscanf(p, "%255s", key) != EOF)
        {
            p += strlen(key);
            while (*p == ' ') p++;
            if (!strcmp(key, "NAME"))
            {
                if (p[strlen(p) - 2] == '\r')
                {
                    p[strlen(p) - 2] = '\0';
                    strcpy(data->name, p);
                }
                if (p[strlen(p) - 1] == '\n')
                {
                    p[strlen(p) - 1] = '\0';
                    strcpy(data->name, p);
                }
                else
                {
                    p[strlen(p) - 2] = '\0';
                    strcpy(data->name, p);
                }
                printf("  Problem Name: %s\n", data->name);
            }
            /*--------------------------------------------------------------------------*/
            /* Type */
            else if (!strcmp(key, "TYPE"))
            {
                data->prob = SOLVER_PROB_OP;
                printf("  Problem Type: %s", p);
                if (sscanf(p, "%255s", field) == EOF)
                {
                    fprintf(stderr, "Missing problem type\n");
                    return 1;
                }
            }
            /*--------------------------------------------------------------------------*/
            /* Comment */
            else if (!strcmp(key, "COMMENT"))
                printf("  Comment: %s", p);
            /*--------------------------------------------------------------------------*/
            /* Dimension */
            else if (!strcmp(key, "DIMENSION"))
            {
                if (sscanf(p, "%255s", field) == EOF)
                {
                    fprintf(stderr, "ERROR in DIMENSION line\n");
                    return 1;
                }
                data->n = atoi(field);
                printf("  Number of Nodes: %d\n", data->n);
            }
        }
    }

    rewind(file);
    while (fgets(buf, 254, file))
    {
        p = buf;
        while (*p != '\0')
        {
            if (*p == ':')
                *p = ' ';
            p++;
        }
        p = buf;
        if (sscanf(p, "%255s", key) != EOF)
        {
            p += strlen(key);
            while (*p == ' ') p++;
            /*--------------------------------------------------------------------------*/
            /* D0 */
            if (!strcmp(key, "COST_LIMIT"))
            {
                if (sscanf(p, "%255s", field) == EOF)
                    printf("Not explicit COST_LIMIT\n");
                data->cap = atof(field);
                printf("  Cost limit: %.2f\n", data->cap);
            }
            else if (!strcmp(key, "DEPOT_SECTION"))
            {
                rval = fscanf(file, "%d", &(data->from));
                check_assert(rval == 1, "", CLEANUP);
                rval = fscanf(file, "%d", &(data->to));
                check_assert(rval == 1, "", CLEANUP);
                data->from--;
                if (data->to == -1)
                    data->to = data->from;
                else
                    data->to--;
                printf("  Depot: %d\n", data->from + 1);
            }
            /*--------------------------------------------------------------------------*/
            /* Edge Weight Type */
            else if (!strcmp(key, "EDGE_WEIGHT_TYPE"))
            {

                if (sscanf(p, "%255s", field) == EOF)
                {
                    fprintf(stderr, "ERROR in EDGE_WEIGHT_TYPE line\n");
                    return 1;
                }
                else if (!strcmp(field, "EXPLICIT"))
                {
                    data->norm = SOLVER_DATA_NORM_MATRIX;
                    printf("  Explicit Lengths (SOLVER_DATA_NORM_MATRIX)\n");
                }
                else if (!strcmp(field, "EUC_2D"))
                {
                    data->norm = SOLVER_DATA_NORM_EUCLIDEAN;
                    printf(
                    "  Rounded Euclidean Norm (SOLVER_DATA_NORM_EUCLIDEAN)\n");
                }
                else if (!strcmp(field, "EUC_3D"))
                {
                    data->norm = SOLVER_DATA_NORM_EUCLIDEAN_3D;
                    printf("  Rounded Euclidean 3D Norm "
                           "(SOLVER_DATA_NORM_EUCLIDEAN_3D)\n");
                }
                else if (!strcmp(field, "GEO"))
                {
                    data->norm = SOLVER_DATA_NORM_GEOGRAPHIC;
                    printf(
                    "  Geographical Norm (SOLVER_DATA_NORM_GEOGRAPHIC)\n");
                }
                else if (!strcmp(field, "GEOM"))
                {
                    data->norm = SOLVER_DATA_NORM_GEOM;
                    printf(
                    "  Geographical Norm in Meters (SOLVER_DATA_NORM_GEOM)\n");
                }
                else if (!strcmp(field, "ATT"))
                {
                    data->norm = SOLVER_DATA_NORM_ATT;
                    printf("  ATT Norm (SOLVER_DATA_NORM_ATT)\n");
                }
                else if (!strcmp(field, "CEIL_2D"))
                {
                    data->norm = SOLVER_DATA_NORM_EUCLIDEAN_CEIL;
                    printf("  Rounded Up Euclidean Norm "
                           "(SOLVER_DATA_NORM_EUCLIDEAN_CEIL)\n");
                }
                else
                {
                    fprintf(stderr, "ERROR: Not set up for norm %s\n", field);
                    return 1;
                }
                data_set_norm_type(data, data->norm);
            }
            /*--------------------------------------------------------------------------*/
            /* Edge Weight Format */
            else if (!strcmp(key, "EDGE_WEIGHT_FORMAT"))
            {
                if (sscanf(p, "%255s", field) == EOF)
                {
                    fprintf(stderr, "ERROR in EDGE_WEIGHT_FORMAT line\n");
                    return 1;
                }
                if (!strcmp(field, "LOWER_DIAG_ROW"))
                    matrixform = MATRIX_LOWER_DIAG_ROW;
                else if (!strcmp(field, "UPPER_ROW"))
                    matrixform = MATRIX_UPPER_ROW;
                else if (!strcmp(field, "UPPER_DIAG_ROW"))
                    matrixform = MATRIX_UPPER_DIAG_ROW;
                else if (!strcmp(field, "FULL_MATRIX"))
                    matrixform = MATRIX_FULL_MATRIX;
                else if (strcmp(field, "FUNCTION"))
                {
                    fprintf(stderr, "Cannot handle format: %s\n", field);
                    return 1;
                }
            }
            /*--------------------------------------------------------------------------*/
            /* Node score format */
            else if (!strcmp(key, "NODE_SCORE_FORMAT"))
            {
                if (sscanf(p, "%255s", field) == EOF)
                {
                    fprintf(stderr, "ERROR in NODE_SCORE_FORMAT line\n");
                    return 1;
                }
                if (!strcmp(field, "FIXED_NODE_SCORES"))
                {
                    printf("Fixed Node scores\n");
                }
                else if (strcmp(field, "FUNCTION"))
                {
                    fprintf(stderr, "Cannot handle format: %s\n", field);
                    return 1;
                }
            }
            /*--------------------------------------------------------------------------*/
            /* Node coordinates section */
            else if (!strcmp(key, "NODE_COORD_SECTION"))
            {
                int i;
                if (data->n <= 0)
                {
                    fprintf(stderr, "ERROR: Dimension not specified\n");
                    return 1;
                }
                if (data->x)
                {
                    fprintf(stderr, "ERROR: A second NODE_COORD_SECTION?\n");
                    data_free(&data);
                    return 1;
                }
                if (data_is_norm_type(data, SOLVER_DATA_TYPE_2D))
                {
                    data->x = malloc(data->n * sizeof(double));
                    if (!data->x)
                    {
                        data_free(&data);
                        return 1;
                    }
                    data->y = malloc(data->n * sizeof(double));
                    if (!data->y)
                    {
                        data_free(&data);
                        return 1;
                    }
                    for (i = 0; i < data->n; i++)
                    {
                        rval = fscanf(file, "%*d %lf %lf", &(data->x[i]),
                                      &(data->y[i]));
                        check_assert(rval == 2, "", CLEANUP);
                    }
                }
                else if (data_is_norm_type(data, SOLVER_DATA_TYPE_3D))
                {
                    data->x = malloc(data->n * sizeof(double));
                    if (!data->x)
                    {
                        data_free(&data);
                        return 1;
                    }
                    data->y = malloc(data->n * sizeof(double));
                    if (!data->y)
                    {
                        data_free(&data);
                        return 1;
                    }
                    data->z = malloc(data->n * sizeof(double));
                    if (!data->z)
                    {
                        data_free(&data);
                        return 1;
                    }
                    for (i = 0; i < data->n; i++)
                    {
                        rval = fscanf(file, "%*d %lf %lf %lf", &(data->x[i]),
                                      &(data->y[i]), &(data->z[i]));
                        check_assert(rval == 3, "", CLEANUP);
                    }
                }
            }
            /*--------------------------------------------------------------------------*/
            /* Node scores section */
            else if (!strcmp(key, "NODE_SCORE_SECTION"))
            {
                int i;
                if (data->n <= 0)
                {
                    fprintf(stderr, "ERROR: Dimension not specified\n");
                    return 1;
                }
                { // op->s = malloc (data->n, sizeof(double));
                    data->obj_node = malloc(data->n * sizeof(double));
                    if (!data->obj_node)
                    {
                        data_free(&data);
                        return 1;
                    }
                    data->tot_obj_node = 0.0;
                    for (i = 0; i < data->n; i++)
                    {
                        rval = fscanf(file, "%*d %lf", &(data->obj_node[i]));
                        check_assert(rval == 1, "", CLEANUP);
                        data->tot_obj_node += data->obj_node[i];
                    }
                }
            }
            /*--------------------------------------------------------------------------*/
            /* Edge weight section */
            else if (!strcmp(key, "EDGE_WEIGHT_SECTION"))
            {
                int i, j;
                if (data->n <= 0)
                {
                    fprintf(stderr, "ERROR: Dimension not specified\n");
                    return 1;
                }
                if (data->adj)
                {
                    fprintf(stderr, "ERROR: A second NODE_COORD_SECTION?\n");
                    data_free(&data);
                    return 1;
                }
                if (data_is_norm_type(data, SOLVER_DATA_TYPE_BANACH))
                {
                    data->adj = malloc(data->n * sizeof(int *));
                    data->adjspace =
                    malloc((data->n) * (data->n + 1) / 2 * sizeof(int));
                    if (!data->adj || !data->adjspace)
                    {
                        data_free(&data);
                        return 1;
                    }
                    for (i = 0, j = 0; i < data->n; i++)
                    {
                        data->adj[i] = data->adjspace + j;
                        j += (i + 1);
                    }
                    if (matrixform == MATRIX_LOWER_DIAG_ROW)
                    {
                        for (i = 0; i < data->n; i++)
                        {
                            for (j = 0; j <= i; j++)
                            {
                                rval = fscanf(file, "%d", &(data->adj[i][j]));
                                check_assert(rval == 1, "", CLEANUP);
                            }
                        }
                    }
                    else if (matrixform == MATRIX_UPPER_ROW ||
                             matrixform == MATRIX_UPPER_DIAG_ROW ||
                             matrixform == MATRIX_FULL_MATRIX)
                    {
                        int **tempadj     = NULL;
                        int *tempadjspace = NULL;
                        tempadj           = malloc(data->n * sizeof(int *));
                        tempadjspace =
                        malloc((data->n) * (data->n) * sizeof(int));
                        if (!tempadj || !tempadjspace)
                        {
                            free(tempadj);
                            free(tempadjspace);
                            data_free(&data);
                            return 1;
                        }
                        for (i = 0; i < data->n; i++)
                        {
                            tempadj[i] = tempadjspace + i * (data->n);
                            if (matrixform == MATRIX_UPPER_ROW)
                            {
                                tempadj[i][i] = 0;
                                for (j = i + 1; j < data->n; j++)
                                {
                                    rval = fscanf(file, "%d", &(tempadj[i][j]));
                                    check_assert(rval == 1, "", CLEANUP);
                                }
                            }
                            else if (matrixform == MATRIX_UPPER_DIAG_ROW)
                            {
                                for (j = i; j < data->n; j++)
                                {
                                    rval = fscanf(file, "%d", &(tempadj[i][j]));
                                    check_assert(rval == 1, "", CLEANUP);
                                }
                            }
                            else
                            {
                                for (j = 0; j < data->n; j++)
                                {
                                    rval = fscanf(file, "%d", &(tempadj[i][j]));
                                    check_assert(rval == 1, "", CLEANUP);
                                }
                            }
                        }
                        for (i = 0; i < data->n; i++)
                        {
                            for (j = 0; j <= i; j++)
                                data->adj[i][j] = tempadj[j][i];
                        }
                        free(tempadjspace);
                        free(tempadj);
                    }
                }
                else
                {
                    fprintf(stderr, "ERROR: Matrix with norm %d?\n",
                            data->norm);
                    return 1;
                }
            }
            /*--------------------------------------------------------------------------*/
            /* Fixed edges section */
            else if (!strcmp(key, "FIXED_EDGES_SECTION"))
            {
                fprintf(stderr, "ERROR: Not set up for fixed edges\n");
                return 1;
            }
        }
    }

CLEANUP:

    if (!data->x && !data->adj)
    {
        fprintf(stderr, "ERROR: Didn't find the data\n");
        return rval;
    }

    if (!rval)
    {
        create_first_map(data);
        data_create_cache(data);
    }
    else
        data_free(&data);

    return rval;
}

static data_map *
create_first_map(solver_data *data)
{
    data_map *map;

    map       = malloc(sizeof(data_map));
    map->fun  = malloc(data->n * sizeof(int));
    map->inv  = malloc(data->n * sizeof(int));
    map->orig = malloc(data->n * sizeof(int));
    if (data->x)
        map->x = malloc(data->n * sizeof(double));
    if (data->y)
        map->y = malloc(data->n * sizeof(double));
    if (data->z)
        map->z = malloc(data->n * sizeof(double));
    else
        map->z = NULL;
    if (data->w)
        map->w = malloc(data->n * sizeof(double));
    else
        map->w = NULL;
    for (int i = 0; i < data->n; i++)
    {
        map->fun[i] = map->inv[i] = map->orig[i] = i;
        if (data->x)
            map->x[i] = data->x[i];
        if (data->y)
            map->y[i] = data->y[i];
        if (data->z)
            map->z[i] = data->z[i];
        if (data->w)
            map->w[i] = data->w[i];
    }
    map->dom_n  = data->n;
    map->img_n  = data->n;
    map->status = 0;
    map->prev   = NULL;
    map->kdtree = NULL;
    data->map   = map;

    if (data_is_norm_type(data, SOLVER_DATA_TYPE_EUCLIDEAN))
        map->kdtree = kdtree_create(data);

    return map;
}

solver_data *
data_read(const char *fname, int format)
{
    FILE *file        = NULL;
    solver_data *data = data_create();

    if (!fname)
    {
        fprintf(stderr, "file name null\n");
        goto done;
    }
    file = fopen(fname, "r");
    if (!file)
    {
        printf("Unable to open '%s'\n", fname);
        goto done;
    }

    printf("Reading problem data from '%s'...\n", fname);

    data_read_(data, file);
    printf("\n");

done:
    if (file)
        fclose(file);
    return data;
}
