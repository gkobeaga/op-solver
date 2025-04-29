#include "op-solver.h"
#include "cp/cp.h"
#include "data/nearest/kdtree/kdtree.h"

struct neighbour
{
    int this;
    int *node;
};

struct position
{
    int prev;
    int next;
    double cost;
};

static void
get_node_3_nearest(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol, int node,
                   struct neighbour *neighbour, int first),
get_best_position(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol, int node,
                  struct neighbour *neighbour, struct position *pos);

int
cp_improve_heur_3n(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol)
{
    int i, prev;
    int nodesel, nodeprev, nodenext;
    int first, completed;
    double len, cost, best, addvalue;
    struct neighbour **neighbours = NULL;
    struct neighbour *neighbour   = NULL;
    struct position *pos          = NULL;
    data_map *map                 = cp->data->map;

    neighbours = malloc(map->img_n * sizeof(struct neighbour *));

    for (i = 0; i < map->img_n; i++)
    {
        if (!sol->selected[i])
        {
            struct neighbour *neighbour = neighbours[i] =
            malloc(sizeof(struct neighbour));
            neighbour->this = i;
            neighbour->node = NULL;
        }
        else
            neighbours[i] = NULL;
    }

    if (map->img_n < 3)
    {
        fprintf(stderr, "Cannot find tour in an %d node graph\n", map->img_n);
        return 1;
    }

    if (data_is_norm_type(cp->data, SOLVER_DATA_TYPE_EUCLIDEAN))
    {
        kdtree_undelete_all(map->kdtree);
        for (i = 0; i < map->img_n; i++)
        {
            if (!sol->selected[i])
                kdtree_delete(map->kdtree, i);
        }
    }

    len       = sol->cap;
    completed = 0;
    first     = 1;

    pos = malloc(sizeof(struct position));

    do
    {
        best = -SOLVER_MAXDOUBLE;
        for (i = 0; i < map->img_n; i++)
        {
            if (!sol->selected[i])
            {
                neighbour = neighbours[i];
                get_node_3_nearest(cp, heur_env, sol, i, neighbour, first);
                get_best_position(cp, heur_env, sol, i, neighbour, pos);

                if (pos->cost <= cp->cap - len)
                    addvalue = cp->data->obj_node[i] / pos->cost;
                else
                    addvalue = -SOLVER_MAXDOUBLE;

                if (addvalue > best)
                {
                    nodesel  = i;
                    nodeprev = pos->prev;
                    nodenext = pos->next;
                    cost     = pos->cost;
                    best     = addvalue;
                }
            }
        }

        if (best > -SOLVER_MAXDOUBLE)
        {
            sol->cod_fr[nodeprev]  = nodesel;
            sol->cod_bk[nodesel]   = nodeprev;
            sol->cod_fr[nodesel]   = nodenext;
            sol->cod_bk[nodenext]  = nodesel;
            sol->selected[nodesel] = 1;
            sol->val += cp->data->obj_node[nodesel];
            sol->ns += 1;

            if (data_is_norm_type(cp->data, SOLVER_DATA_TYPE_EUCLIDEAN))
                kdtree_undelete(map->kdtree, nodesel);

            for (i = 0; i < map->img_n; i++)
            {
                if (!sol->selected[i])
                {
                    neighbour          = neighbours[i];
                    neighbour->node[3] = nodesel;
                }
            }

            len += cost;
            first     = 0;
            neighbour = neighbours[nodesel];
            free(neighbour->node);
            neighbour->node = NULL;
        }
        else
        {
            completed = 1;
        }

    } while (completed == 0 && sol->ns < map->img_n);

    prev = cp->data->from;
    for (i = 0; i < sol->ns; i++)
    {
        sol->cycle[i] = prev;
        prev          = sol->cod_fr[prev];
    }
    assert(prev == cp->data->from);
    for (i = sol->ns; i < map->img_n; i++) sol->cycle[i] = -1;

    sol->cap = len;

    if (data_is_norm_type(cp->data, SOLVER_DATA_TYPE_EUCLIDEAN))
        kdtree_undelete_all(map->kdtree);

    for (i = 0; i < map->img_n; i++)
    {
        neighbour = neighbours[i];
        if (neighbour)
        {
            free(neighbour->node);
            free(neighbour);
        }
    }
    free(neighbours);
    free(pos);
    return 0;
}

int
cp_add_sol_node(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol, int node)
{
    int rval = 0;
    int i, j, prev;
    struct neighbour *neighbour = NULL;
    struct position *pos        = NULL;
    data_map *map               = cp->data->map;

    neighbour       = malloc(sizeof(struct neighbour));
    neighbour->this = node;
    neighbour->node = malloc(map->img_n * sizeof(int));

    if (map->img_n < 3)
    {
        fprintf(stderr, "Cannot find tour in an %d node graph\n", map->img_n);
    }

    if (data_is_norm_type(cp->data, SOLVER_DATA_TYPE_EUCLIDEAN))
    {
        kdtree_undelete_all(map->kdtree);
        for (i = 1; i < map->img_n; i++)
        {
            if (!sol->selected[i])
                kdtree_delete(map->kdtree, i);
        }
    }

    pos = malloc(sizeof(struct position));

    get_node_3_nearest(cp, heur_env, sol, node, neighbour, 1);
    get_best_position(cp, heur_env, sol, node, neighbour, pos);

    sol->cod_fr[pos->prev] = node;
    sol->cod_bk[node]      = pos->prev;
    sol->cod_fr[node]      = pos->next;
    sol->cod_bk[pos->next] = node;
    sol->cap += pos->cost;
    sol->val += cp->data->obj_node[node];
    sol->selected[node] = 1;
    sol->ns += 1;

    prev = 0;
    for (i = 0; i < sol->ns; i++)
    {
        sol->cycle[i] = prev;
        prev          = sol->cod_fr[prev];
    }
    assert(prev == 0);
    for (i = sol->ns; i < map->img_n; i++) sol->cycle[i] = -1;

    j = 0;
    for (i = 0; i < map->img_n; i++)
    {
        if (sol->selected[i])
        {
            sol->sposition[j] = i;
            j++;
        }
    }

    if (data_is_norm_type(cp->data, SOLVER_DATA_TYPE_EUCLIDEAN))
        kdtree_undelete_all(map->kdtree);

    free(neighbour->node);
    free(neighbour);
    free(pos);
    return rval;
}

static void
get_node_3_nearest(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol, int node,
                   struct neighbour *neighbour, int first)
{
    int i, k;
    data_map *map = cp->data->map;

    if (first || sol->ns < 4)
    {
        if (sol->ns > 3 &&
            data_is_norm_type(cp->data, SOLVER_DATA_TYPE_EUCLIDEAN))
        {
            data_get_node_k_nearest(cp->data, node, 3, &neighbour->node);
            neighbour->node = realloc(neighbour->node, 4 * sizeof(int));
        }
        else
        {
            neighbour->node = malloc(4 * sizeof(int));
            double *len     = malloc(4 * sizeof(double));
            double edgelen;

            for (i = 0; i < 3; i++) len[i] = SOLVER_MAXDOUBLE;
            len[3] = 0.0;

            for (i = 0; i < map->img_n; i++)
            {
                if (sol->selected[i])
                {
                    edgelen = (double)data_get_norm(cp->data, node, i);
                    if (len[0] > edgelen)
                    {
                        for (k = 0; len[k + 1] > edgelen; k++)
                        {
                            neighbour->node[k] = neighbour->node[k + 1];
                            len[k]             = len[k + 1];
                        }
                        neighbour->node[k] = i;
                        len[k]             = edgelen;
                    }
                }
            }

            free(len);
        }
    }
    else
    {
        double *len = malloc(4 * sizeof(double));
        double edgelen;

        for (i = 0; i < 3; i++) len[i] = SOLVER_MAXDOUBLE;
        len[3] = 0.0;

        for (i = 0; i < map->img_n; i++)
        {
            if (node != i &&
                (i == neighbour->node[0] || i == neighbour->node[1] ||
                 i == neighbour->node[2] || i == neighbour->node[3]))
            {
                edgelen = (double)data_get_norm(cp->data, node, i);
                if (len[0] > edgelen)
                {
                    for (k = 0; len[k + 1] > edgelen; k++)
                    {
                        neighbour->node[k] = neighbour->node[k + 1];
                        len[k]             = len[k + 1];
                    }
                    neighbour->node[k] = i;
                    len[k]             = edgelen;
                }
            }
        }

        free(len);
    }
}

void
get_best_position(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol, int node,
                  struct neighbour *neighbour, struct position *pos)
{
    int prev, next;
    double cost, best = SOLVER_MAXDOUBLE;

    if (sol->ns > 2)
    {
        if (sol->cod_fr[neighbour->node[0]] == neighbour->node[1] ||
            sol->cod_fr[neighbour->node[1]] == neighbour->node[0] ||
            sol->cod_fr[neighbour->node[0]] == neighbour->node[2] ||
            sol->cod_fr[neighbour->node[2]] == neighbour->node[0] ||
            sol->cod_fr[neighbour->node[1]] == neighbour->node[2] ||
            sol->cod_fr[neighbour->node[2]] == neighbour->node[1])
        {
            if (sol->cod_fr[neighbour->node[0]] == neighbour->node[1])
            {
                best =
                (double)(data_get_norm(cp->data, neighbour->node[0], node) +
                         data_get_norm(cp->data, node, neighbour->node[1]) -
                         data_get_norm(cp->data, neighbour->node[0],
                                       neighbour->node[1]));
                prev = neighbour->node[0];
                next = neighbour->node[1];
            }
            if (sol->cod_fr[neighbour->node[1]] == neighbour->node[0])
            {
                cost =
                (double)(data_get_norm(cp->data, neighbour->node[1], node) +
                         data_get_norm(cp->data, node, neighbour->node[0]) -
                         data_get_norm(cp->data, neighbour->node[1],
                                       neighbour->node[0]));
                if (cost < best)
                {
                    prev = neighbour->node[1];
                    next = neighbour->node[0];
                    best = cost;
                }
            }
            if (sol->cod_fr[neighbour->node[0]] == neighbour->node[2])
            {
                cost =
                (double)(data_get_norm(cp->data, neighbour->node[0], node) +
                         data_get_norm(cp->data, node, neighbour->node[2]) -
                         data_get_norm(cp->data, neighbour->node[0],
                                       neighbour->node[2]));
                if (cost < best)
                {
                    prev = neighbour->node[0];
                    next = neighbour->node[2];
                    best = cost;
                }
            }
            if (sol->cod_fr[neighbour->node[2]] == neighbour->node[0])
            {
                cost =
                (double)(data_get_norm(cp->data, neighbour->node[2], node) +
                         data_get_norm(cp->data, node, neighbour->node[0]) -
                         data_get_norm(cp->data, neighbour->node[2],
                                       neighbour->node[0]));
                if (cost < best)
                {
                    prev = neighbour->node[2];
                    next = neighbour->node[0];
                    best = cost;
                }
            }
            if (sol->cod_fr[neighbour->node[1]] == neighbour->node[2])
            {
                cost =
                (double)(data_get_norm(cp->data, neighbour->node[1], node) +
                         data_get_norm(cp->data, node, neighbour->node[2]) -
                         data_get_norm(cp->data, neighbour->node[1],
                                       neighbour->node[2]));
                if (cost < best)
                {
                    prev = neighbour->node[1];
                    next = neighbour->node[2];
                    best = cost;
                }
            }
            if (sol->cod_fr[neighbour->node[2]] == neighbour->node[1])
            {
                cost =
                (double)(data_get_norm(cp->data, neighbour->node[2], node) +
                         data_get_norm(cp->data, node, neighbour->node[1]) -
                         data_get_norm(cp->data, neighbour->node[2],
                                       neighbour->node[1]));
                if (cost < best)
                {
                    prev = neighbour->node[2];
                    next = neighbour->node[1];
                    best = cost;
                }
            }
        }
        else
        {
            best = (double)(data_get_norm(cp->data, neighbour->node[0], node) +
                            data_get_norm(cp->data, node,
                                          sol->cod_fr[neighbour->node[0]]) -
                            data_get_norm(cp->data, neighbour->node[0],
                                          sol->cod_fr[neighbour->node[0]]));
            prev = neighbour->node[0];
            next = sol->cod_fr[neighbour->node[0]];

            cost = (double)(data_get_norm(cp->data, neighbour->node[1], node) +
                            data_get_norm(cp->data, node,
                                          sol->cod_fr[neighbour->node[1]]) -
                            data_get_norm(cp->data, neighbour->node[1],
                                          sol->cod_fr[neighbour->node[1]]));
            if (cost < best)
            {
                prev = neighbour->node[1];
                next = sol->cod_fr[neighbour->node[1]];
                best = cost;
            }

            cost = (double)(data_get_norm(cp->data, neighbour->node[2], node) +
                            data_get_norm(cp->data, node,
                                          sol->cod_fr[neighbour->node[2]]) -
                            data_get_norm(cp->data, neighbour->node[2],
                                          sol->cod_fr[neighbour->node[2]]));
            if (cost < best)
            {
                prev = neighbour->node[2];
                next = sol->cod_fr[neighbour->node[2]];
                best = cost;
            }
        }
        pos->prev = prev;
        pos->next = next;
        pos->cost = best;
    }
    else if (sol->ns == 1)
    {
        pos->prev = 0;
        pos->next = 0;
        pos->cost = (double)(2 * data_get_norm(cp->data, 0, node));
    }
    else if (sol->ns == 2)
    {
        pos->prev = 0;
        pos->next = sol->cod_fr[0];
        pos->cost = (double)(data_get_norm(cp->data, 0, node) +
                             data_get_norm(cp->data, node, sol->cod_fr[0]) -
                             data_get_norm(cp->data, 0, sol->cod_fr[0]));
    }
}

int
cp_improve_heur_in(cp_prob *cp, cp_heur_env *heur_env, cp_sol *sol)
{
    int rval = 0;
    int i, j, node, prev, next;
    int try = 1;
    double min, cap_tmp, value, tmp_value;
    data_map *map = cp->data->map;

    while (try)
    {
        value = -SOLVER_MAXDOUBLE;
        for (i = 0; i < map->img_n; i++)
        {
            if (!sol->selected[i])
            {
                for (j = 0; j < sol->ns; j++)
                {
                    if (j < sol->ns - 1)
                    {
                        cap_tmp =
                        (double)(data_get_norm(cp->data, sol->cycle[j], i) +
                                 data_get_norm(cp->data, i, sol->cycle[j + 1]) -
                                 data_get_norm(cp->data, sol->cycle[j],
                                               sol->cycle[j + 1]));
                    }
                    else
                    {
                        cap_tmp =
                        (double)(data_get_norm(cp->data, sol->cycle[j], i) +
                                 data_get_norm(cp->data, i, sol->cycle[0]) -
                                 data_get_norm(cp->data, sol->cycle[j],
                                               sol->cycle[0]));
                    }
                    if (sol->cap + cap_tmp < cp->cap)
                        tmp_value = cp->data->obj_node[i] / cap_tmp;
                    else
                        tmp_value = -SOLVER_MAXDOUBLE;

                    if (tmp_value > value)
                    {
                        min   = cap_tmp;
                        value = tmp_value;
                        prev  = sol->cycle[j];
                        if (j < sol->ns - 1)
                            next = sol->cycle[j + 1];
                        else
                            next = sol->cycle[0];
                        node = i;
                    }
                }
            }
        }

        if (value > -SOLVER_MAXDOUBLE)
        {
            sol->cod_fr[prev] = node;
            sol->cod_bk[node] = prev;
            sol->cod_fr[node] = next;
            sol->cod_bk[next] = node;
            sol->cap += min;
            sol->val += cp->data->obj_node[node];
            sol->selected[node] = 1;
            sol->ns += 1;

            prev = 0;
            for (i = 0; i < sol->ns; i++)
            {
                sol->cycle[i] = prev;
                prev          = sol->cod_fr[prev];
            }
            check_rval(prev == 0, "ERROR: prev!=0\n", CLEANUP);
        }
        else
            try = 0;
    }

    for (i = 0, j = 0; i < map->img_n; i++)
    {
        if (sol->selected[i])
            sol->sposition[j++] = i;
    }

CLEANUP:

    return rval;
}
