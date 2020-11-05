#include "op-solver.h"

unsigned int
clique_hash(graph_clique *c)
{
    unsigned int x = 0;
    int i;

    for (i = 0; i < c->segcount; i++)
    {
        graph_clique_segment seg = c->nodes[i];

        x = x * 65537 + seg.lo * 4099 + seg.hi;
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
