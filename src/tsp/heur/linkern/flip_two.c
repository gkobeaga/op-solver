#include "op-solver.h"
#include "tsp/tsp.h"
#include "linkern.h"

#define GROUPSIZE_FACTOR 0.50
#define SEGMENT_SPLIT_CUTOFF 0.30

static void
same_segment_flip(__lk_flipper *F, __lk_childnode *a, __lk_childnode *b),
consecutive_segment_flip(__lk_flipper *F, __lk_parentnode *a,
                         __lk_parentnode *b),
segment_split(__lk_flipper *F, __lk_parentnode *p, __lk_childnode *aprev,
              __lk_childnode *a, int left_or_right),
init_flipper(__lk_flipper *Fl), free_flipper(__lk_flipper *Fl);

static int
build_flipper(__lk_flipper *Fl, int ncount);

#define SAME_SEGMENT(a, b)                                                     \
    (a->parent == b->parent &&                                                 \
     ((!((F->reversed) ^ (a->parent->rev)) && a->id <= b->id) ||               \
      (((F->reversed) ^ (a->parent->rev)) && a->id >= b->id)))

int
__linkern_flipper_init(__lk_flipper *F, tsp_sol *sol)
{
    int i, j, cind, remain;
    int rval = 0;
    __lk_childnode *c, *cprev;
    __lk_parentnode *p;

    init_flipper(F);
    rval = build_flipper(F, sol->tot_n);
    if (rval)
    {
        fprintf(stderr, "build_flipper failed\n");
        goto CLEANUP;
    }

    remain = sol->tot_n;
    i      = 0;
    j      = 2 * F->groupsize;
    while (remain >= j)
    {
        F->parents[i].size = F->groupsize;
        remain -= F->groupsize;
        i++;
    }
    if (remain > F->groupsize)
    {
        F->parents[i].size = remain / 2;
        remain -= (remain / 2);
        i++;
    }
    F->parents[i].size = remain;
    i++;

    if (i != F->nsegments)
    {
        fprintf(stderr, "seg count is wrong\n");
        rval = 1;
        goto CLEANUP;
    }

    c = &(F->children[sol->cycle[sol->tot_n - 1]]);
    for (i = 0, p = F->parents, cind = 0; i < F->nsegments; p++, i++)
    {
        p->id      = i;
        p->rev     = 0;
        p->ends[0] = &(F->children[sol->cycle[cind]]);
        for (j = p->size; j > 0; j--)
        {
            cprev         = c;
            c             = &(F->children[sol->cycle[cind]]);
            c->id         = cind;
            c->name       = sol->cycle[cind];
            c->parent     = p;
            c->adj[0]     = cprev;
            cprev->adj[1] = c;
            cind++;
        }
        p->ends[1] = c;
        p->adj[0]  = p - 1;
        p->adj[1]  = p + 1;
    }
    F->parents[0].adj[0]                = &(F->parents[F->nsegments - 1]);
    F->parents[F->nsegments - 1].adj[1] = &(F->parents[0]);

CLEANUP:

    if (rval)
    {
        free_flipper(F);
    }
    return rval;
}

void
__linkern_flipper_cycle(__lk_flipper *F, tsp_sol *sol, solver_data *data)
{
    __lk_childnode *c, *start;
    int k = 0, prev;

    start = &(F->children[0]);
    c     = start->adj[!((F->reversed) ^ (start->parent->rev))];

    tsp_erase_sol(sol);

    sol->cycle[k++]            = start->name;
    sol->selected[start->name] = 1;
    sol->ns++;
    prev = start->name;
    while (c != start)
    {
        sol->cycle[k++]        = c->name;
        sol->cod_fr[prev]      = c->name;
        sol->cod_bk[c->name]   = prev;
        sol->selected[c->name] = 1;
        sol->ns++;
        sol->val += data_get_norm(data, prev, c->name);
        prev = c->name;
        c    = c->adj[!((F->reversed) ^ (c->parent->rev))];
    }
    sol->cod_fr[prev]    = c->name;
    sol->cod_bk[c->name] = prev;
    sol->val += data_get_norm(data, prev, c->name);
}

void
__linkern_flipper_finish(__lk_flipper *F)
{
    free_flipper(F);
}

int
__linkern_flipper_next(__lk_flipper *F, int x)
{
    return F->children[x]
    .adj[!((F->reversed) ^ (F->children[x].parent->rev))]
    ->name;
}

int
__linkern_flipper_prev(__lk_flipper *F, int x)
{
    return F->children[x]
    .adj[(F->reversed) ^ (F->children[x].parent->rev)]
    ->name;
}

void
__linkern_flipper_flip(__lk_flipper *F, int x, int y)
{
    __lk_childnode *xc = &(F->children[x]);
    __lk_childnode *yc = &(F->children[y]);

    if (SAME_SEGMENT(xc, yc))
    {
        if (xc != yc)
        {
            same_segment_flip(F, xc, yc);
        }
    }
    else
    {
        int xdir              = ((F->reversed) ^ (xc->parent->rev));
        int ydir              = ((F->reversed) ^ (yc->parent->rev));
        __lk_childnode *xprev = xc->adj[xdir];
        __lk_childnode *ynext = yc->adj[!ydir];
        if (SAME_SEGMENT(ynext, xprev))
        {
            if (ynext != xprev)
            {
                same_segment_flip(F, ynext, xprev);
            }
            (F->reversed) ^= 1;
        }
        else
        {
            int side;
            if (xc->parent->ends[xdir] == xc && yc->parent->ends[!ydir] == yc)
            {
                if (F->reversed)
                    side = xc->parent->id - yc->parent->id;
                else
                    side = yc->parent->id - xc->parent->id;
                if (side < 0)
                    side += F->nsegments;
                if (side < F->nsegments / 2)
                {
                    consecutive_segment_flip(F, xc->parent, yc->parent);
                }
                else
                {
                    consecutive_segment_flip(F, yc->parent->adj[!F->reversed],
                                             xc->parent->adj[F->reversed]);
                    (F->reversed) ^= 1;
                }
            }
            else
            {
                if (xprev->parent == xc->parent)
                {
                    segment_split(F, xc->parent, xprev, xc, 0);
                    if (SAME_SEGMENT(xc, yc))
                    {
                        if (xc != yc)
                            same_segment_flip(F, xc, yc);
                        return;
                    }
                    else if (SAME_SEGMENT(ynext, xprev))
                    {
                        if (ynext != xprev)
                        {
                            same_segment_flip(F, ynext, xprev);
                        }
                        (F->reversed) ^= 1;
                        return;
                    }
                }
                if (ynext->parent == yc->parent)
                {
                    segment_split(F, yc->parent, yc, ynext, 0);
                    if (SAME_SEGMENT(xc, yc))
                    {
                        if (xc != yc)
                            same_segment_flip(F, xc, yc);
                        return;
                    }
                    else if (SAME_SEGMENT(ynext, xprev))
                    {
                        if (ynext != xprev)
                        {
                            same_segment_flip(F, ynext, xprev);
                        }
                        (F->reversed) ^= 1;
                        return;
                    }
                }
                if (F->reversed)
                    side = xc->parent->id - yc->parent->id;
                else
                    side = yc->parent->id - xc->parent->id;
                if (side < 0)
                    side += F->nsegments;
                if (side < F->nsegments / 2)
                {
                    consecutive_segment_flip(F, xc->parent, yc->parent);
                }
                else
                {
                    consecutive_segment_flip(F, yc->parent->adj[!F->reversed],
                                             xc->parent->adj[F->reversed]);
                    (F->reversed) ^= 1;
                }
            }
        }
    }
}

static void
same_segment_flip(__lk_flipper *F, __lk_childnode *a, __lk_childnode *b)
{
    __lk_parentnode *parent = a->parent;
    int dir                 = ((F->reversed) ^ (parent->rev));
    __lk_childnode *aprev   = a->adj[dir];
    __lk_childnode *bnext   = b->adj[!dir];
    __lk_childnode *c, *cnext;

    if ((dir && a->id - b->id > F->split_cutoff) ||
        (!dir && b->id - a->id > F->split_cutoff))
    {
        if (aprev->parent == parent)
            segment_split(F, parent, aprev, a, 1);
        if (bnext->parent == parent)
            segment_split(F, parent, b, bnext, 2);
        aprev->adj[!((F->reversed) ^ (aprev->parent->rev))] = b;
        bnext->adj[(F->reversed) ^ (bnext->parent->rev)]    = a;
        a->adj[dir]                                         = bnext;
        b->adj[!dir]                                        = aprev;
        parent->rev ^= 1;
        return;
    }

    if (dir)
    {
        int id                                              = a->id;
        aprev->adj[!((F->reversed) ^ (aprev->parent->rev))] = b;
        bnext->adj[(F->reversed) ^ (bnext->parent->rev)]    = a;
        cnext                                               = b->adj[1];
        b->adj[1]                                           = aprev;
        b->adj[0]                                           = cnext;
        b->id                                               = id--;
        c                                                   = cnext;
        while (c != a)
        {
            cnext     = c->adj[1];
            c->adj[1] = c->adj[0];
            c->adj[0] = cnext;
            c->id     = id--;
            c         = cnext;
        }
        a->adj[1] = a->adj[0];
        a->adj[0] = bnext;
        a->id     = id;
        if (parent->ends[1] == a)
            parent->ends[1] = b;
        if (parent->ends[0] == b)
            parent->ends[0] = a;
    }
    else
    {
        int id                                              = a->id;
        aprev->adj[!((F->reversed) ^ (aprev->parent->rev))] = b;
        bnext->adj[(F->reversed) ^ (bnext->parent->rev)]    = a;
        c                                                   = b->adj[0];
        b->adj[0]                                           = aprev;
        b->adj[1]                                           = c;
        b->id                                               = id++;
        while (c != a)
        {
            cnext     = c->adj[0];
            c->adj[0] = c->adj[1];
            c->adj[1] = cnext;
            c->id     = id++;
            c         = cnext;
        }
        a->adj[0] = a->adj[1];
        a->adj[1] = bnext;
        a->id     = id;
        if (parent->ends[0] == a)
            parent->ends[0] = b;
        if (parent->ends[1] == b)
            parent->ends[1] = a;
    }
}

static void
consecutive_segment_flip(__lk_flipper *F, __lk_parentnode *a,
                         __lk_parentnode *b)
{
    __lk_parentnode *aprev = a->adj[F->reversed];
    __lk_parentnode *bnext = b->adj[!F->reversed];
    __lk_parentnode *c, *cnext;
    __lk_childnode *achild = a->ends[(F->reversed) ^ (a->rev)];
    __lk_childnode *bchild = b->ends[!((F->reversed) ^ (b->rev))];
    __lk_childnode *childprev, *childnext;
    int id = a->id;

    if (F->reversed)
    {
        childprev                               = achild->adj[!a->rev];
        childnext                               = bchild->adj[b->rev];
        childprev->adj[childprev->parent->rev]  = bchild;
        childnext->adj[!childnext->parent->rev] = achild;
        bchild->adj[b->rev]                     = childprev;
        achild->adj[!a->rev]                    = childnext;

        aprev->adj[0] = b;
        bnext->adj[1] = a;
        c             = b->adj[1];
        b->adj[1]     = aprev;
        b->adj[0]     = c;
        b->id         = id--;
        b->rev ^= 1;
        while (c != a)
        {
            cnext     = c->adj[1];
            c->adj[1] = c->adj[0];
            c->adj[0] = cnext;
            c->id     = id--;
            c->rev ^= 1;
            c = cnext;
        }
        a->adj[1] = a->adj[0];
        a->adj[0] = bnext;
        a->id     = id;
        a->rev ^= 1;
    }
    else
    {
        childprev                               = achild->adj[a->rev];
        childnext                               = bchild->adj[!b->rev];
        childprev->adj[!childprev->parent->rev] = bchild;
        childnext->adj[childnext->parent->rev]  = achild;
        bchild->adj[!b->rev]                    = childprev;
        achild->adj[a->rev]                     = childnext;

        aprev->adj[1] = b;
        bnext->adj[0] = a;
        c             = b->adj[0];
        b->adj[0]     = aprev;
        b->adj[1]     = c;
        b->id         = id++;
        b->rev ^= 1;
        while (c != a)
        {
            cnext     = c->adj[0];
            c->adj[0] = c->adj[1];
            c->adj[1] = cnext;
            c->id     = id++;
            c->rev ^= 1;
            c = cnext;
        }
        a->adj[0] = a->adj[1];
        a->adj[1] = bnext;
        a->id     = id;
        a->rev ^= 1;
    }
}

/* split between a and aprev */
static void
segment_split(__lk_flipper *F, __lk_parentnode *p, __lk_childnode *aprev,
              __lk_childnode *a, int left_or_right)
{
    int side;
    int dir = ((F->reversed) ^ (p->rev));
    int id;
    __lk_parentnode *pnext;
    __lk_childnode *b, *bnext;

    if (dir)
        side = p->ends[1]->id - aprev->id + 1;
    else
        side = aprev->id - p->ends[0]->id + 1;

    if ((left_or_right == 0 && side <= p->size / 2) || left_or_right == 1)
    {
        pnext = p->adj[F->reversed];
        pnext->size += side;
        p->size -= side;
        if (pnext->rev == p->rev)
        {
            b  = pnext->ends[!dir];
            id = b->id;
            if (dir)
            {
                do
                {
                    b         = b->adj[0];
                    b->id     = --id;
                    b->parent = pnext;
                } while (b != aprev);
            }
            else
            {
                do
                {
                    b         = b->adj[1];
                    b->id     = ++id;
                    b->parent = pnext;
                } while (b != aprev);
            }
            pnext->ends[!dir] = aprev;
            p->ends[dir]      = a;
        }
        else
        {
            b  = pnext->ends[dir];
            id = b->id;
            if (!dir)
            {
                bnext = b->adj[0];
                do
                {
                    b         = bnext;
                    b->id     = --id;
                    b->parent = pnext;
                    bnext     = b->adj[1];
                    b->adj[1] = b->adj[0];
                    b->adj[0] = bnext;
                } while (b != aprev);
            }
            else
            {
                bnext = b->adj[1];
                do
                {
                    b         = bnext;
                    b->id     = ++id;
                    b->parent = pnext;
                    bnext     = b->adj[0];
                    b->adj[0] = b->adj[1];
                    b->adj[1] = bnext;
                } while (b != aprev);
            }
            pnext->ends[dir] = aprev;
            p->ends[dir]     = a;
        }
    }
    else
    {
        pnext = p->adj[!F->reversed];
        pnext->size += (p->size - side);
        p->size = side;
        if (pnext->rev == p->rev)
        {
            b  = pnext->ends[dir];
            id = b->id;
            if (dir)
            {
                do
                {
                    b         = b->adj[1];
                    b->id     = ++id;
                    b->parent = pnext;
                } while (b != a);
            }
            else
            {
                do
                {
                    b         = b->adj[0];
                    b->id     = --id;
                    b->parent = pnext;
                } while (b != a);
            }
            pnext->ends[dir] = a;
            p->ends[!dir]    = aprev;
        }
        else
        {
            b  = pnext->ends[!dir];
            id = b->id;
            if (!dir)
            {
                bnext = b->adj[1];
                do
                {
                    b         = bnext;
                    b->id     = ++id;
                    b->parent = pnext;
                    bnext     = b->adj[0];
                    b->adj[0] = b->adj[1];
                    b->adj[1] = bnext;
                } while (b != a);
            }
            else
            {
                bnext = b->adj[0];
                do
                {
                    b         = bnext;
                    b->id     = --id;
                    b->parent = pnext;
                    bnext     = b->adj[1];
                    b->adj[1] = b->adj[0];
                    b->adj[0] = bnext;
                } while (b != a);
            }
            pnext->ends[!dir] = a;
            p->ends[!dir]     = aprev;
        }
    }
}

int
__linkern_flipper_sequence(__lk_flipper *F, int x, int y, int z)
{
    __lk_childnode *a   = &(F->children[x]);
    __lk_childnode *b   = &(F->children[y]);
    __lk_childnode *c   = &(F->children[z]);
    __lk_parentnode *pa = a->parent;
    __lk_parentnode *pb = b->parent;
    __lk_parentnode *pc = c->parent;

    if (pa == pb)
    {
        if (pa == pc)
        {
            if ((F->reversed) ^ (pa->rev))
            {
                if (a->id >= b->id)
                {
                    return (b->id >= c->id || c->id >= a->id);
                }
                else
                {
                    return (b->id >= c->id && c->id >= a->id);
                }
            }
            else
            {
                if (a->id <= b->id)
                {
                    return (b->id <= c->id || c->id <= a->id);
                }
                else
                {
                    return (b->id <= c->id && c->id <= a->id);
                }
            }
        }
        else
        {
            if ((F->reversed) ^ (pa->rev))
            {
                return (a->id >= b->id);
            }
            else
            {
                return (a->id <= b->id);
            }
        }
    }
    else if (pa == pc)
    {
        if ((F->reversed) ^ (pa->rev))
        {
            return (a->id <= c->id);
        }
        else
        {
            return (a->id >= c->id);
        }
    }
    else if (pb == pc)
    {
        if ((F->reversed) ^ (pb->rev))
        {
            return (b->id >= c->id);
        }
        else
        {
            return (b->id <= c->id);
        }
    }
    else
    {
        if (F->reversed)
        {
            if (pa->id >= pb->id)
            {
                return (pb->id >= pc->id || pc->id >= pa->id);
            }
            else
            {
                return (pb->id >= pc->id && pc->id >= pa->id);
            }
        }
        else
        {
            if (pa->id <= pb->id)
            {
                return (pb->id <= pc->id || pc->id <= pa->id);
            }
            else
            {
                return (pb->id <= pc->id && pc->id <= pa->id);
            }
        }
    }
}

static void
init_flipper(__lk_flipper *Fl)
{
    Fl->parents      = NULL;
    Fl->children     = NULL;
    Fl->reversed     = 0;
    Fl->nsegments    = 0;
    Fl->groupsize    = 100;
    Fl->split_cutoff = 100;
}

static void
free_flipper(__lk_flipper *Fl)
{
    if (Fl)
    {
        free(Fl->parents);
        free(Fl->children);
        Fl->reversed     = 0;
        Fl->nsegments    = 0;
        Fl->groupsize    = 0;
        Fl->split_cutoff = 0;
    }
}

static int
build_flipper(__lk_flipper *Fl, int ncount)
{
    int rval = 0;

    Fl->reversed     = 0;
    Fl->groupsize    = (int)(sqrt((double)ncount) * GROUPSIZE_FACTOR);
    Fl->nsegments    = (ncount + Fl->groupsize - 1) / Fl->groupsize;
    Fl->split_cutoff = Fl->groupsize * SEGMENT_SPLIT_CUTOFF;

    Fl->parents  = malloc(Fl->nsegments * sizeof(__lk_parentnode));
    Fl->children = malloc((ncount + 1) * sizeof(__lk_childnode));
    /* The +1 will stop a purify burp later */
    if (!Fl->parents || !Fl->children)
    {
        fprintf(stderr, "out of memory in build_flipper\n");
        rval = 1;
        goto CLEANUP;
    }

CLEANUP:

    if (rval)
    {
        free_flipper(Fl);
    }
    return rval;
}
