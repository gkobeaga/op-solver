#include "op-solver.h"

graph_clique *
clique_create(void)
{
    graph_clique *clique = malloc(sizeof(graph_clique));
    clique->segcount     = 0;
    clique->nodes        = NULL;
    clique->hashnext     = 0;
    clique->refcount     = 0;
    clique->val          = 0.0;
    clique->prev         = NULL;
    clique->next         = NULL;
    return clique;
}

static int
sort_vert(const void *xx, const void *yy)
{
    graph_vertex *x = *(graph_vertex **)xx, *y = *(graph_vertex **)yy;

    if (x->i < y->i)
        return -1;
    if (x->i > y->i)
        return +1;

    return 0;
}

graph_clique *
clique_conv_vertices2clique(solver_graph *graph, graph_vertex **vert, int vcount)
{
    graph_clique *clique = clique_create();
    int i, nseg;

    qsort(vert, vcount, sizeof(graph_vertex *), sort_vert);

    nseg = 0;

    i = 0;
    while (i < vcount)
    {
        assert(i == (vcount - 1) || vert[i + 1]->i != vert[i]->i);
        while (i < (vcount - 1) && vert[i + 1]->i == (vert[i]->i + 1)) i++;
        i++;
        nseg++;
    }

    clique->nodes = malloc(nseg * sizeof(graph_clique_segment));
    if (!clique->nodes)
    {
        fprintf(stderr, "out of memory in conv_vertices2clique\n");
        free(clique);
        return NULL;
    }
    clique->segcount = nseg;

    nseg = 0;
    i    = 0;
    while (i < vcount)
    {
        clique->nodes[nseg].lo = vert[i]->i;
        while (i < (vcount - 1) && vert[i + 1]->i == (vert[i]->i + 1)) i++;
        clique->nodes[nseg].hi = vert[i++]->i;
        nseg++;
    }

    return clique;
}

graph_clique *
clique_conv_vertices2coclique(solver_graph *graph, graph_vertex **vert, int vcount)
{
    graph_clique *clique = clique_create();
    int i, nseg;

    qsort(vert, vcount, sizeof(graph_vertex *), sort_vert);

    nseg = 0;

    i = 0;
    while (i < vcount)
    {
        assert(i == (vcount - 1) || vert[i + 1]->i != vert[i]->i);
        if (0 < vert[i]->i)
            nseg++;
        while (i < (vcount - 1) && vert[i + 1]->i == (vert[i]->i + 1)) i++;
        i++;
    }
    if (vert[vcount - 1]->i < graph->nv - 1)
        nseg++;

    clique->nodes = malloc(nseg * sizeof(graph_clique_segment));
    if (!clique->nodes)
    {
        fprintf(stderr, "out of memory in conv_vertices2coclique\n");
        free(clique);
        return NULL;
    }
    clique->segcount = nseg;

    nseg = 0;
    i    = 0;
    if (vert[i]->i)
    {
        clique->nodes[nseg].lo = 0;
        clique->nodes[nseg].hi = vert[i]->i - 1;
        nseg++;
    }
    while (i < (vcount - 1) && vert[i + 1]->i == (vert[i]->i + 1)) i++;
    while (i < vcount - 1)
    {
        clique->nodes[nseg].lo = vert[i]->i + 1;
        clique->nodes[nseg].hi = vert[++i]->i - 1;
        nseg++;
        while (i < (vcount - 1) && vert[i + 1]->i == (vert[i]->i + 1)) i++;
    }
    if (vert[i]->i < graph->nv - 1)
    {
        clique->nodes[nseg].lo = vert[i]->i + 1;
        clique->nodes[nseg].hi = graph->nv - 1;
        nseg++;
    }

    return clique;
}

static int
sort_array(const void *xx, const void *yy)
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

int
clique_conv_array2clique(int *vert, int vcount, graph_clique **clique)
{
    int rval = 0;
    *clique  = clique_create();
    int i, nseg;

    qsort(vert, vcount, sizeof(int), sort_array);

    nseg = 0;
    i    = 0;
    while (i < vcount)
    {
        while (i < (vcount - 1) && vert[i + 1] == (vert[i] + 1)) i++;
        i++;
        nseg++;
    }

    (*clique)->nodes = malloc(nseg * sizeof(graph_clique_segment));
    check_null((*clique)->nodes, "out of memory", CLEANUP);
    (*clique)->segcount = nseg;

    nseg = 0;
    i    = 0;
    while (i < vcount)
    {
        (*clique)->nodes[nseg].lo = vert[i];
        while (i < (vcount - 1) && vert[i + 1] == (vert[i] + 1)) i++;
        (*clique)->nodes[nseg].hi = vert[i++];
        nseg++;
    }

CLEANUP:

    return rval;
}

int
clique_conv_clique2array(graph_clique *clique, int **ar, int *count)
{
    int rval = 0;
    int j, tmp;
    int k = 0;

    *ar = (int *)NULL;

    *count = clique_count(clique);
    if (*count)
    {
        *ar = malloc(*count * sizeof(int));
        check_null(*ar, "out of memory", CLEANUP);

        FOREACH_NODE_IN_CLIQUE (j, clique, tmp)
        {
            (*ar)[k++] = j;
        }
    }

CLEANUP:

    return rval;
}

int
clique_copy(graph_clique *in, graph_clique *out)
{
    int rval = 0;
    int k;
    graph_clique_segment *s = NULL;

    if (in->segcount)
    {
        s = malloc(in->segcount * sizeof(graph_clique_segment));
        check_null(s, "out of memory", CLEANUP);

        for (k = 0; k < in->segcount; k++)
        {
            s[k].lo = in->nodes[k].lo;
            s[k].hi = in->nodes[k].hi;
        }
    }
    out->segcount = in->segcount;
    out->val      = in->val;
    out->nodes    = s;
CLEANUP:

    return rval;
}

void
clique_free(graph_clique **clique)
{
    if (*clique)
    {
        if ((*clique)->nodes)
        {
            free((*clique)->nodes);
        }
        free(*clique);
        *clique = NULL;
    }
}

unsigned int
clique_hash(graph_clique *c)
{
    unsigned int x = 0;
    int i;

    for (i = 0; i < c->segcount; i++)
    {
        graph_clique_segment seg = c->nodes[i];
        x                        = x * 65537 + seg.lo * 4099 + seg.hi;
    }

    return x;
}

int
clique_eq(graph_clique *c, graph_clique *d)
{
    int i;

    if (c->segcount != d->segcount)
        return 0;
    for (i = 0; i < c->segcount; i++)
    {
        graph_clique_segment seg1 = c->nodes[i];
        graph_clique_segment seg2 = d->nodes[i];
        if (seg1.lo != seg2.lo)
            return 0;
        if (seg1.hi != seg2.hi)
            return 0;
    }
    return 1;
}

int
clique_count(graph_clique *clique)
{
    int i;
    int count                   = 0;
    graph_clique_segment *nodes = clique->nodes;

    for (i = 0; i < clique->segcount; i++)
        count += (nodes[i].hi - nodes[i].lo + 1);
    return count;
}

void
clique_print(graph_clique *clique)
{
    int i, tmp;
    int count1, count2 = 0;

    printf("Clique: (count %d)\n", count1 = clique_count(clique));
    FOREACH_NODE_IN_CLIQUE (i, clique, tmp)
    {
        count2++;
        printf(" %d", i);
    }
    printf("\n");
    assert(count1 == count2);

    return;
}
