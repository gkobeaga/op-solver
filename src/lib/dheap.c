#include "op-solver.h"

#define dheap_D 3
#define dheap_UP(x) (((x) - 1) / dheap_D)
#define dheap_DOWN(x) (((x) * dheap_D) + 1)

static void
dheap_siftup(solver_dheap *dheap, int i, int x),
dheap_siftdown(solver_dheap *dheap, int i, int x);

static int
dheap_minchild(solver_dheap *dheap, int x);

solver_dheap *
dheap_create(int k)
{
    solver_dheap *dheap = malloc(sizeof(solver_dheap));
    dheap->loc          = NULL;
    dheap->key          = NULL;
    dheap->entry        = malloc(k * sizeof(int));
    if (!dheap->entry)
        return NULL;
    dheap->loc = malloc(k * sizeof(int));
    if (!dheap->loc)
    {
        free(dheap->entry);
        return NULL;
    }
    dheap->key = malloc(k * sizeof(double));
    if (!dheap->key)
    {
        free(dheap->entry);
        free(dheap->loc);
        return NULL;
    }
    dheap->total_space = k;
    dheap->size        = 0;
    return dheap;
}

void
dheap_free(solver_dheap **dheap)
{
    if (*dheap)
    {
        free((*dheap)->entry);
        free((*dheap)->loc);
        free((*dheap)->key);
        free(*dheap);
        *dheap = NULL;
    }
}

int
dheap_resize(solver_dheap *dheap, int newsize)
{
    int rval = 0;
    if (newsize < dheap->size || newsize < dheap->total_space)
        return rval;
    dheap->key = realloc(dheap->key, newsize);
    check_null(dheap->key, "trealloc failed", CLEANUP);
    dheap->entry = realloc(dheap->entry, newsize);
    check_null(dheap->entry, "trealloc failed", CLEANUP);
    dheap->loc = realloc(dheap->loc, newsize);
    check_null(dheap->loc, "trealloc failed", CLEANUP);
    dheap->total_space = newsize;

    return rval;

CLEANUP:
    return rval;
}

int
dheap_findmin(solver_dheap *dheap)
{
    if (dheap->size == 0)
        return -1;
    else
        return dheap->entry[0];
}

int
dheap_insert(solver_dheap *dheap, int i)
{
    if (dheap->size >= dheap->total_space)
    {
        fprintf(stderr, "Error - dheap already full\n");
        return 1;
    }
    dheap->size++;
    dheap_siftup(dheap, i, dheap->size - 1);
    return 0;
}

void
dheap_delete(solver_dheap *dheap, int i)
{
    int j;

    dheap->size--;
    j                         = dheap->entry[dheap->size];
    dheap->entry[dheap->size] = -1;

    if (j != i)
    {
        if (dheap->key[j] <= dheap->key[i])
            dheap_siftup(dheap, j, dheap->loc[i]);
        else
            dheap_siftdown(dheap, j, dheap->loc[i]);
    }
}

int
dheap_deletemin(solver_dheap *dheap)
{
    int i;

    if (dheap->size == 0)
        return -1;
    else
    {
        i = dheap->entry[0];
        dheap_delete(dheap, i);
        return i;
    }
}

void
dheap_changekey(solver_dheap *dheap, int i, double newkey)
{
    if (newkey < dheap->key[i])
    {
        dheap->key[i] = newkey;
        dheap_siftup(dheap, i, dheap->loc[i]);
    }
    else if (newkey > dheap->key[i])
    {
        dheap->key[i] = newkey;
        dheap_siftdown(dheap, i, dheap->loc[i]);
    }
}

static void
dheap_siftup(solver_dheap *dheap, int i, int x)
{
    int p;

    p = dheap_UP(x);
    while (x && dheap->key[dheap->entry[p]] > dheap->key[i])
    {
        dheap->entry[x]             = dheap->entry[p];
        dheap->loc[dheap->entry[p]] = x;
        x                           = p;
        p                           = dheap_UP(p);
    }
    dheap->entry[x] = i;
    dheap->loc[i]   = x;
}

static void
dheap_siftdown(solver_dheap *dheap, int i, int x)
{
    int c;

    c = dheap_minchild(dheap, x);

    while (c >= 0 && dheap->key[dheap->entry[c]] < dheap->key[i])
    {
        dheap->entry[x]             = dheap->entry[c];
        dheap->loc[dheap->entry[c]] = x;
        x                           = c;
        c                           = dheap_minchild(dheap, c);
    }
    dheap->entry[x] = i;
    dheap->loc[i]   = x;
}

static int
dheap_minchild(solver_dheap *dheap, int x)
{
    int c = dheap_DOWN(x);
    int cend;
    double minval;
    int minloc;

    if (c >= dheap->size)
        return -1;
    minval = dheap->key[dheap->entry[c]];
    minloc = c;
    cend   = c + dheap_D;
    if (dheap->size < cend)
        cend = dheap->size;
    for (c++; c < cend; c++)
    {
        if (dheap->key[dheap->entry[c]] < minval)
        {
            minval = dheap->key[dheap->entry[c]];
            minloc = c;
        }
    }
    return minloc;
}
