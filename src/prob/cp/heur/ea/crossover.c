#include "op-solver.h"
#include "cp/heur/ea/ea.h"

#define BIGINT 2000000000

struct selnode
{
    int nextsel1;
    int nextsel2;
    int prevsel1;
    int prevsel2;
    int *next_internodes1;
    int *next_internodes2;
    int *prev_internodes1;
    int *prev_internodes2;
    int next_ninter1;
    int next_ninter2;
    int prev_ninter1;
    int prev_ninter2;
};

static void
initselnodes(struct selnode *selnode),
freeselnodes(struct selnode *selnode),
reset_relative_degree(int ncount, int *relative_degree, double *unvisited,
                      int *degree, struct selnode *selnodes, int current),
update_neigh_degrees(int *degree, struct selnode *selnodes, int current),
insert_inter_nodes(struct selnode *selnodes, int current, int next,
                   cp_sol *child);

static int
cp_heur_ea_crossover_work(cp_prob *cp, cp_heur_ea_env *ea_env, cp_sol *par0,
                          cp_sol *par1, cp_sol *child),
select_connected_next(int ncount, struct selnode *selnodes, int current,
                      double *unvisited, int *relative_degree),
select_disconnected_next(int ncount, int current, double *unvisited,
                         int *degree),
sort_int(const void *xx, const void *yy),
all_degree_zero(int ncount, int *degree);

int
cp_heur_ea_crossover(cp_prob *cp, cp_heur_ea_env *ea_env, cp_sol *par0,
                     cp_sol *par1, cp_sol *child)
{
    int rval = 0;
    int i, ncommon;

    ncommon = 0;
    for (i = 0; i < cp->n; i++)
    {
        if (par0->cod_fr[i] == par1->cod_fr[i])
            ncommon++;
    }

    if (ncommon != par0->ns || ncommon != par1->ns)
    {
        rval = cp_heur_ea_crossover_work(cp, ea_env, par0, par1, child);
        check_rval(rval, "crossover failed", cleanup);

        child->val = 0.0;
        for (i = 0; i < cp->n; i++)
        {
            if (child->selected[i])
                child->val += cp->data->obj_node[i];
        }
        child->cap = (double)data_get_norm(
        cp->data, child->cycle[child->ns - 1], child->cycle[0]);
        for (i = 1; i < child->ns; i++)
            child->cap += (double)data_get_norm(cp->data, child->cycle[i - 1],
                                                child->cycle[i]);
    }
    else
    {
        cp_sol *par;
        if (rand() % 2)
            par = par0;
        else
            par = par1;
        cp_copy_sol(par, child);
    }

cleanup:
    return rval;
}

static int
cp_heur_ea_crossover_work(cp_prob *cp, cp_heur_ea_env *ea_env, cp_sol *par0,
                          cp_sol *par1, cp_sol *child)
{
    int rval = 0;
    int i, j;
    int next1, next2;
    int *relative_degree = malloc(cp->n * sizeof(int));
    int node;
    int bscount, nunvisited;
    int *tinternodes  = NULL;
    int *degree_init  = NULL;
    int *degree       = NULL;
    double *unvisited = NULL;
    int *neightnodes  = NULL;
    int ninter1, ninter2;
    int current, next;
    struct selnode *selnodes;

    // We select the nodes that are in the two paths
    bscount = 0; // count nodes visited in both paths
    for (i = 0; i < cp->n; i++)
    {
        child->selected[i] = par0->selected[i] * par1->selected[i];
        if (child->selected[i])
            bscount++;
    }

    selnodes = malloc(cp->n * sizeof(struct selnode));

    for (i = 0; i < cp->n; i++) initselnodes(&selnodes[i]);

    // Temp array to save intermediate nodes
    tinternodes = malloc(cp->n * sizeof(int));
    check_null(tinternodes, "out of memory", cleanup);

    for (i = 0; i < cp->n; i++)
    {
        if (child->selected[i])
        { // For parent 1
            ninter1 = 0;
            next1   = par0->cod_fr[i];
            while (!child->selected[next1])
            {
                tinternodes[ninter1] = next1;
                ninter1 += 1;

                next1 = par0->cod_fr[next1];
            }

            selnodes[i].nextsel1         = next1;
            selnodes[i].next_ninter1     = ninter1;
            selnodes[next1].prevsel1     = i;
            selnodes[next1].prev_ninter1 = ninter1;

            if (ninter1 != 0)
            {
                selnodes[i].next_internodes1 = malloc(ninter1 * sizeof(int));
                selnodes[next1].prev_internodes1 =
                malloc(ninter1 * sizeof(int));
                for (j = 0; j < ninter1; j++)
                {
                    selnodes[i].next_internodes1[j] = tinternodes[j];
                    selnodes[next1].prev_internodes1[ninter1 - j - 1] =
                    tinternodes[j];
                }
            }

            // For parent 2
            ninter2 = 0;
            next2   = par1->cod_fr[i];
            while (!child->selected[next2])
            {
                tinternodes[ninter2] = next2;
                ninter2 += 1;

                next2 = par1->cod_fr[next2];
            }

            selnodes[i].nextsel2         = next2;
            selnodes[i].next_ninter2     = ninter2;
            selnodes[next2].prevsel2     = i;
            selnodes[next2].prev_ninter2 = ninter2;

            if (ninter2 != 0)
            {
                selnodes[i].next_internodes2 = malloc(ninter2 * sizeof(int));
                selnodes[next2].prev_internodes2 =
                malloc(ninter2 * sizeof(int));
                for (j = 0; j < ninter2; j++)
                {
                    selnodes[i].next_internodes2[j] = tinternodes[j];
                    selnodes[next2].prev_internodes2[ninter2 - j - 1] =
                    tinternodes[j];
                }
            }
        }
    }

    // Save the vertex degree to use in the crossover. How many diferent
    // neighbours have in the the reduced parents.

    degree_init = malloc(cp->n * sizeof(int));
    check_null(degree_init, "out of memory", cleanup);

    for (i = 0; i < cp->n; i++)
    {
        if (child->selected[i])
        {
            neightnodes = malloc(4 * sizeof(int));
            check_null(neightnodes, "out of memory", cleanup);

            neightnodes[0] = selnodes[i].nextsel1;
            neightnodes[1] = selnodes[i].prevsel1;
            neightnodes[2] = selnodes[i].nextsel2;
            neightnodes[3] = selnodes[i].prevsel2;

            qsort(neightnodes, 4, sizeof(int), sort_int);

            degree_init[i] = 1;
            node           = neightnodes[0];
            for (j = 1; j < 4; j++)
            {
                if (neightnodes[j] != node)
                {
                    degree_init[i] += 1;
                    node = neightnodes[j];
                }
            }

            free(neightnodes);
        }
        else
        {
            degree_init[i] = BIGINT;
        }
    }

    degree = malloc(cp->n * sizeof(int));
    check_null(degree, "out of memory", cleanup);
    memset(degree, 0, cp->n * sizeof(int));
    for (i = 0; i < cp->n; i++)
    {
        if (child->selected[i])
            degree[i] = degree_init[i];
    }

    nunvisited = bscount;
    unvisited  = malloc(cp->n * sizeof(double));
    check_null(unvisited, "out of memory", cleanup);
    memset(unvisited, 0.0, cp->n * sizeof(double));
    for (i = 1; i < cp->n; i++)
    {
        if (child->selected[i])
            unvisited[i] = 1.0;
    }

    child->ns               = 0;
    current                 = 0;
    child->cycle[child->ns] = current;
    child->ns++;

    unvisited[current] = 0.0;
    nunvisited--;
    update_neigh_degrees(degree, selnodes, current);

    while (nunvisited > 0)
    {
        reset_relative_degree(cp->n, relative_degree, unvisited, degree,
                              selnodes, current);

        if (all_degree_zero(cp->n, relative_degree))
        {
            next = select_connected_next(cp->n, selnodes, current, unvisited,
                                         relative_degree);

            insert_inter_nodes(selnodes, current, next, child);

            current                 = next;
            child->cycle[child->ns] = current;
            child->ns++;
            unvisited[current] = 0.0;
            nunvisited--;
            update_neigh_degrees(degree, selnodes, current);

            // If hasn't got unvisited connected nodes. Select one unvisited
            // node randomly as next, and another connected (visited) node for
            // inserting intermediates nodes.
        }
        else
        {
            next = select_disconnected_next(cp->n, current, unvisited, degree);

            child->cod_fr[current] = next;
            child->cod_bk[next]    = current;

            current                 = next;
            child->cycle[child->ns] = current;
            child->ns++;
            unvisited[current] = 0;
            nunvisited--;
            update_neigh_degrees(degree, selnodes, current);
        }
    }

    // if current and 0 connected
    if ((selnodes[current].nextsel1 == 0 || selnodes[current].prevsel1 == 0) &&
        (selnodes[current].nextsel2 == 0 || selnodes[current].prevsel2 == 0))
    {
        insert_inter_nodes(selnodes, current, 0, child);
    }
    else
    {
        child->cod_fr[current] = 0;
        child->cod_bk[0]       = current;
    }

    j = 0;
    for (i = 0; i < cp->n; i++)
    {
        if (child->selected[i])
        {
            child->sposition[j] = i;
            j++;
        }
        else
        {
            child->cod_bk[i] = i;
            child->cod_fr[i] = i;
        }
    }

    for (i = child->ns; i < cp->n; i++) child->cycle[i] = -1;

cleanup:
    free(degree);
    free(relative_degree);
    free(unvisited);
    free(degree_init);
    free(tinternodes);
    for (i = 0; i < cp->n; i++) freeselnodes(&selnodes[i]);
    free(selnodes);
    return rval;
}

static int
sort_int(const void *xx, const void *yy)
/**************************************************************************/
{
    int x = *((int *)xx);
    int y = *((int *)yy);

    if (x < y)
        return -1;
    if (x > y)
        return +1;

    return 0;
}

static int
all_degree_zero(int ncount, int *degree)
/**************************************************************************/
{
    int i;

    for (i = 0; i < ncount; i++)
    {
        if (degree[i] != 0)
            return 1;
    }

    return 0;
}

static void
reset_relative_degree(int ncount, int *relative_degree, double *unvisited,
                      int *degree, struct selnode *selnodes, int current)
{
    memset(relative_degree, 0, ncount * sizeof(int));
    if (unvisited[selnodes[current].nextsel1])
        relative_degree[selnodes[current].nextsel1] =
        degree[selnodes[current].nextsel1];
    if (unvisited[selnodes[current].prevsel1])
        relative_degree[selnodes[current].prevsel1] =
        degree[selnodes[current].prevsel1];
    if (unvisited[selnodes[current].nextsel2])
        relative_degree[selnodes[current].nextsel2] =
        degree[selnodes[current].nextsel2];
    if (unvisited[selnodes[current].prevsel2])
        relative_degree[selnodes[current].prevsel2] =
        degree[selnodes[current].prevsel2];
}

static void
update_neigh_degrees(int *degree, struct selnode *selnodes, int current)
{
    int i, node;
    int *neightnodes = malloc(4 * sizeof(int));

    neightnodes[0] = selnodes[current].nextsel1;
    neightnodes[1] = selnodes[current].prevsel1;
    neightnodes[2] = selnodes[current].nextsel2;
    neightnodes[3] = selnodes[current].prevsel2;

    qsort(neightnodes, 4, sizeof(int), sort_int);

    node = neightnodes[0];
    degree[node]--;
    for (i = 1; i < 4; i++)
    {
        if (neightnodes[i] != node)
        {
            node = neightnodes[i];
            degree[node]--;
        }
    }

    free(neightnodes);
}

static int
select_connected_next(int ncount, struct selnode *selnodes, int current,
                      double *unvisited, int *relative_degree)
{
    int i, j, next;
    double min, val;
    int *candidates;

    min = BIGINT;

    for (i = 0; i < ncount; i++)
    {
        val = relative_degree[i];
        if (val != 0 && val < min)
            min = val;
    }

    // Step 3 of the generalized edge recombination algorithm
    if (min != BIGINT)
    {
        candidates = malloc(ncount * sizeof(int));

        for (i = 0, j = 0; i < ncount; i++)
        {
            if (relative_degree[i] == min)
                candidates[j++] = i;
        }

        next = candidates[rand() % j];

        free(candidates);

        return next;
    }
    else
    {
        // If current hasn't got a unvisited connected neighbour. We return 0.
        // Notice that 0 never could be the next node, because we start the
        // algorithm in node 0.
        // TODO: update this if a alternative starting node is considered
        return 0;
    }
}

static int
select_disconnected_next(int ncount, int current, double *unvisited,
                         int *degree)
{
    int i, j, next;
    int *candidates;

    candidates = malloc(ncount * sizeof(int));

    for (i = 0, j = 0; i < ncount; i++)
    {
        if (unvisited[i])
            candidates[j++] = i;
    }

    next = candidates[rand() % j];

    free(candidates);

    return next;
}

// This function might be written in another way using the information
// in cycle1 and cycle2.
static void
insert_inter_nodes(struct selnode *selnodes, int current, int next,
                   cp_sol *child)
{
    int i;
    int parent, node;

    // If there are multiple connection between current and next nodes.
    if ((selnodes[current].nextsel1 == next ||
         selnodes[current].prevsel1 == next) &&
        (selnodes[current].nextsel2 == next ||
         selnodes[current].prevsel2 == next))
    {
        node   = current;
        parent = rand() % 2;
        if (parent == 0)
        {
            if (selnodes[current].nextsel1 == next &&
                selnodes[current].next_ninter1 != 0)
            {
                for (i = 0; i < selnodes[current].next_ninter1; i++)
                {
                    child->cod_fr[node] = selnodes[current].next_internodes1[i];
                    child->cod_bk[selnodes[current].next_internodes1[i]] = node;
                    node = selnodes[current].next_internodes1[i];
                    child->cycle[child->ns] = node;
                    child->selected[node]   = 1;
                    child->ns++;
                }
            }
            else if (selnodes[current].prevsel1 == next &&
                     selnodes[current].prev_ninter1 != 0)
            {
                for (i = 0; i < selnodes[current].prev_ninter1; i++)
                {
                    child->cod_fr[node] = selnodes[current].prev_internodes1[i];
                    child->cod_bk[selnodes[current].prev_internodes1[i]] = node;
                    node = selnodes[current].prev_internodes1[i];
                    child->cycle[child->ns] = node;
                    child->selected[node]   = 1;
                    child->ns++;
                }
            }
        }
        else
        {
            node = current;
            if (selnodes[current].nextsel2 == next &&
                selnodes[current].next_ninter2 != 0)
            {
                for (i = 0; i < selnodes[current].next_ninter2; i++)
                {
                    child->cod_fr[node] = selnodes[current].next_internodes2[i];
                    child->cod_bk[selnodes[current].next_internodes2[i]] = node;
                    node = selnodes[current].next_internodes2[i];
                    child->cycle[child->ns] = node;
                    child->selected[node]   = 1;
                    child->ns++;
                }
            }
            else if (selnodes[current].prevsel2 == next &&
                     selnodes[current].prev_ninter2 != 0)
            {
                for (i = 0; i < selnodes[current].prev_ninter2; i++)
                {
                    child->cod_fr[node] = selnodes[current].prev_internodes2[i];
                    child->cod_bk[selnodes[current].prev_internodes2[i]] = node;
                    node = selnodes[current].prev_internodes2[i];
                    child->cycle[child->ns] = node;
                    child->selected[node]   = 1;
                    child->ns++;
                }
            }
        }
        child->cod_fr[node] = next;
        child->cod_bk[next] = node;
        // If there is a unique connection between current and next nodes.
    }
    else
    {
        node = current;
        if (selnodes[current].nextsel1 == next &&
            selnodes[current].next_ninter1 != 0)
        {
            for (i = 0; i < selnodes[current].next_ninter1; i++)
            {
                child->cod_fr[node] = selnodes[current].next_internodes1[i];
                child->cod_bk[selnodes[current].next_internodes1[i]] = node;
                node                    = selnodes[current].next_internodes1[i];
                child->cycle[child->ns] = node;
                child->selected[node]   = 1;
                child->ns++;
            }
        }
        else if (selnodes[current].prevsel1 == next &&
                 selnodes[current].prev_ninter1 != 0)
        {
            for (i = 0; i < selnodes[current].prev_ninter1; i++)
            {
                child->cod_fr[node] = selnodes[current].prev_internodes1[i];
                child->cod_bk[selnodes[current].prev_internodes1[i]] = node;
                node                    = selnodes[current].prev_internodes1[i];
                child->cycle[child->ns] = node;
                child->selected[node]   = 1;
                child->ns++;
            }
        }
        else if (selnodes[current].nextsel2 == next &&
                 selnodes[current].next_ninter2 != 0)
        {
            for (i = 0; i < selnodes[current].next_ninter2; i++)
            {
                child->cod_fr[node] = selnodes[current].next_internodes2[i];
                child->cod_bk[selnodes[current].next_internodes2[i]] = node;
                node                    = selnodes[current].next_internodes2[i];
                child->cycle[child->ns] = node;
                child->selected[node]   = 1;
                child->ns += 1;
            }
        }
        else if (selnodes[current].prevsel2 == next &&
                 selnodes[current].prev_ninter2 != 0)
        {
            for (i = 0; i < selnodes[current].prev_ninter2; i++)
            {
                child->cod_fr[node] = selnodes[current].prev_internodes2[i];
                child->cod_bk[selnodes[current].prev_internodes2[i]] = node;
                node                    = selnodes[current].prev_internodes2[i];
                child->cycle[child->ns] = node;
                child->selected[node]   = 1;
                child->ns += 1;
            }
        }
        child->cod_fr[node] = next;
        child->cod_bk[next] = node;
    }
}

static void
initselnodes(struct selnode *selnode)
{
    selnode->nextsel1         = 0;
    selnode->nextsel2         = 0;
    selnode->prevsel1         = 0;
    selnode->prevsel2         = 0;
    selnode->next_internodes1 = (int *)NULL;
    selnode->next_internodes2 = (int *)NULL;
    selnode->prev_internodes1 = (int *)NULL;
    selnode->prev_internodes2 = (int *)NULL;
    selnode->next_ninter1     = 0;
    selnode->next_ninter2     = 0;
    selnode->prev_ninter1     = 0;
    selnode->prev_ninter2     = 0;
}

static void
freeselnodes(struct selnode *selnode)
{
    if (selnode)
    {
        selnode->nextsel1 = 0;
        selnode->nextsel2 = 0;
        selnode->prevsel1 = 0;
        selnode->prevsel2 = 0;
        if (selnode->next_internodes1 != (int *)NULL)
            free(selnode->next_internodes1);
        if (selnode->next_internodes2 != (int *)NULL)
            free(selnode->next_internodes2);
        if (selnode->prev_internodes1 != (int *)NULL)
            free(selnode->prev_internodes1);
        if (selnode->prev_internodes2 != (int *)NULL)
            free(selnode->prev_internodes2);
        selnode->next_ninter1 = 0;
        selnode->next_ninter2 = 0;
        selnode->prev_ninter1 = 0;
        selnode->prev_ninter2 = 0;
    }
}
