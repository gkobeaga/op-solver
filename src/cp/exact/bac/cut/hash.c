#include "cp/cp.h"
#include "ip/exact/bac/bac.h"
#include "cp/exact/bac/bac.h"

static int
genhash_resize(cp_cut_hash *h);

#define CUTHASH_MAXDENSITY 1.0
#define CUTHASH_LOWDENSITY 0.6

cp_cut_repo *
cp_create_cut_repo(cp_prob *cp)
{
    cp_cut_repo *repo = malloc(sizeof(cp_cut_repo));
    repo->count       = 0;
    repo->space       = 0;
    repo->cliques     = clique_create_repo(cp->n);
    repo->hash        = cp_create_cut_hash(
    (int)((SOLVER_IP_BAC_PRICE_GEN + SOLVER_IP_BAC_PRICE_GEN_FACTOR * cp->n) *
          1.5));
    repo->cuts = NULL;
    return repo;
}

void
cp_free_cut_repo(cp_cut_repo **repo)
{
    if (*repo)
    {
        free((*repo)->hash->table);
        free((*repo)->hash);
        clique_free_repo(&(*repo)->cliques);
        if ((*repo)->cuts)
        {
            for (int i = 0; i < (*repo)->count; i++)
                cp_free_cut(&(*repo)->cuts[i]);
            free((*repo)->cuts);
        }
        free((*repo));
        *repo = NULL;
    }
}

int
create_cut_hash_work(cp_cut_hash *hash, int size)
{
    int i;

    hash->nelem      = 0;
    hash->size       = prime_next((unsigned int)size);
    hash->maxdensity = CUTHASH_MAXDENSITY;
    hash->lowdensity = CUTHASH_LOWDENSITY;

    hash->maxelem = (int)(hash->maxdensity * hash->size);

    hash->table = malloc(hash->size * sizeof(cp_cut *));
    if (!hash->table)
        return -1;

    for (i = 0; i < hash->size; i++) hash->table[i] = NULL;

    return 0;
}

cp_cut_hash *
cp_create_cut_hash(int size)
{
    cp_cut_hash *hash = malloc(sizeof(cp_cut_hash));

    create_cut_hash_work(hash, size);
    return hash;
}

unsigned int
cp_get_cut_hash(cp_cut *cut)
{
    unsigned int x =
    ((unsigned int)cut->rhs) * 257 + ((unsigned int)cut->sense);
    int i, j;

    if (cut->hcount + cut->tcount)
    {
        for (i = 0; i < cut->hcount; i++) x = x * 4099 + cut->handle_cid[i];
        for (i = 0; i < cut->tcount; i++) x = x * 4099 + cut->teeth_cid[i];

        for (i = 0; i < 2 * cut->tcount; i++)
        {
            if (cut->verts[i] >= 0)
                x = x * 4099 + cut->verts[i];
        }
    }
    else if (cut->logical)
    {
        x = x * 4099 + cut->logical->arc[0];
        x = x * 4099 + cut->logical->arc[1];
        x = x * 4099 + cut->logical->v;
    }
    else if (cut->cover_edge)
    {
        for (i = 0; i < 2 * cut->cover_edge->na; i++)
            x = x * 4099 + cut->cover_edge->arcs[i];
    }
    else if (cut->connect)
    {
        FOREACH_NODE_IN_CLIQUE (i, cut->connect->verts, j)
        {
            x = x * 4099 + i;
        }
    }
    else if (cut->path)
    {
        for (i = 0; i < 2 * cut->path->na; i++)
            x = x * 4099 + cut->path->arcs[i];
        for (i = 0; i < 2 * cut->path->fna; i++)
            x = x * 4099 + cut->path->farcs[i];
    }

    return x;
}

int
cp_eq_cut(cp_cut *cut1, cp_cut *cut2)
{
    int i;

    if (cut1->hcount + cut1->tcount && cut2->hcount + cut2->tcount)
    {
        if (cut1->hcount != cut2->hcount)
            return 1;
        if (cut1->tcount != cut2->tcount)
            return 1;
        if (cut1->rhs != cut2->rhs)
            return 1;
        if (cut1->sense != cut2->sense)
            return 1;
        for (i = 0; i < 2 * cut1->tcount; i++)
            if (cut1->verts[i] != cut2->verts[i])
                return 1;
        if (cut1->vycoef != cut2->vycoef)
            return 1;
        for (i = 0; i < cut1->hcount; i++)
        {
            if (cut1->handle_cid[i] != cut2->handle_cid[i])
                return 1;
        }
        for (i = 0; i < cut1->tcount; i++)
        {
            if (cut1->teeth_cid[i] != cut2->teeth_cid[i])
                return 1;
        }
    }
    else if (cut1->logical && cut2->logical)
    {
        if (cut1->logical->arc[0] != cut2->logical->arc[0])
            return 1;
        if (cut1->logical->arc[1] != cut2->logical->arc[1])
            return 1;
        if (cut1->logical->v != cut2->logical->v)
            return 1;
    }
    else if (cut1->cover_edge && cut2->cover_edge)
    {
        if (cut1->cover_edge->na != cut2->cover_edge->na)
            return 1;
        if (cut1->cover_edge->strong != cut2->cover_edge->strong)
            return 1;
        for (i = 0; i < 2 * cut1->cover_edge->na; i++)
            if (cut1->cover_edge->arcs[i] != cut2->cover_edge->arcs[i])
                return 1;
    }
    else if (cut1->connect && cut2->connect)
    {
        if (cut1->connect->verts->segcount != cut2->connect->verts->segcount)
            return 1;
        for (i = 0; i < cut1->connect->verts->segcount; i++)
        {
            if (cut1->connect->verts->nodes[i].lo !=
                cut2->connect->verts->nodes[i].lo)
                return 1;
            if (cut1->connect->verts->nodes[i].hi !=
                cut2->connect->verts->nodes[i].hi)
                return 1;
        }
    }
    else if (cut1->path && cut2->path)
    {
        if (cut1->path->na != cut2->path->na)
            return 1;
        if (cut1->path->fna != cut2->path->fna)
            return 1;
        for (i = 0; i < 2 * cut1->path->na; i++)
            if (cut1->path->arcs[i] != cut2->path->arcs[i])
                return 1;
        for (i = 0; i < 2 * cut1->path->fna; i++)
            if (cut1->path->farcs[i] != cut2->path->farcs[i])
                return 1;
    }
    else
        return 1;

    return 0;
}

int
cp_add_cut_hash(cp_cut_hash *hash, cp_cut *cut)
{
    int loc;
    unsigned int hashval;

    if (hash->nelem >= hash->maxelem)
    {
        if (genhash_resize(hash))
        {
            return -1;
        }
    }

    hashval = cp_get_cut_hash(cut);
    loc     = hashval % hash->size;

    cut->hash_next = hash->table[loc];
    if (cut->hash_next)
        cut->hash_next->hash_prev = cut;
    hash->table[loc] = cut;
    return 0;
}

int
cp_del_cut_hash(cp_cut_hash *hash, cp_cut *cut)
{
    int loc = cp_get_cut_hash(cut) % hash->size;
    cp_cut *tmpcut;

    for (tmpcut = hash->table[loc]; tmpcut; tmpcut = tmpcut->hash_next)
    {
        if (!cp_eq_cut(cut, tmpcut))
        {
            if (tmpcut->hash_prev)
                tmpcut->hash_prev->hash_next = tmpcut->hash_next;
            else
                hash->table[loc] = tmpcut->hash_next;
            if (tmpcut->hash_next)
                tmpcut->hash_next->hash_prev = tmpcut->hash_prev;
            return 0;
        }
    }
    return 1;
}

cp_cut *
cp_find_cut_hash(cp_cut_hash *hash, cp_cut *cut)
{
    cp_cut *tmpcut;
    int loc = cp_get_cut_hash(cut) % hash->size;

    for (tmpcut = hash->table[loc]; tmpcut; tmpcut = tmpcut->hash_next)
    {
        if (!cp_eq_cut(cut, tmpcut))
        {
            return tmpcut;
        }
    }
    return NULL;
}

static int
genhash_resize(cp_cut_hash *h)
{
    int newsize = prime_next((unsigned int)(h->nelem / h->lowdensity));
    cp_cut **newtable;
    cp_cut *cut, *nextcut;
    int loc;
    int i;

    if (newsize <= h->nelem)
        newsize = prime_next((unsigned int)h->nelem + 1);

    newtable = malloc(newsize * sizeof(cp_cut **));
    if (!newtable)
        return -1;

    for (i = 0; i < newsize; i++) newtable[i] = NULL;

    for (i = 0; i < h->size; i++)
    {
        for (cut = h->table[i]; cut; cut = nextcut)
        {
            nextcut       = cut->hash_next;
            loc           = cp_get_cut_hash(cut) % newsize;
            cut->next     = newtable[loc];
            newtable[loc] = cut;
        }
    }

    free(h->table);

    h->table   = newtable;
    h->size    = newsize;
    h->maxelem = (int)(h->maxdensity * h->size);

    return 0;
}

static int
cp_register_cut_repo_clique(cp_cut_repo *cuts, graph_clique *c)
{
    int x = clique_hash(c) % cuts->cliques->cliquehashsize;
    int y = cuts->cliques->cliquehash[x];
    int i;
    graph_clique *clique;

    while (y != -1)
    {
        clique = cuts->cliques->cliques[y];
        if (clique_eq(c, clique))
        {
            clique->refcount++;
            return y;
        }
        y = clique->hashnext;
    }

    if (cuts->cliques->cliquefree != -1)
    {
        y                         = cuts->cliques->cliquefree;
        clique                    = cuts->cliques->cliques[y];
        cuts->cliques->cliquefree = clique->hashnext;
        clique                    = cuts->cliques->cliques[y];
    }
    else
    {
        if (cuts->cliques->size + 1 >= cuts->cliques->cliquespace)
        {
            realloc_scale(cuts->cliques->cliques, cuts->cliques->cliquespace,
                          cuts->cliques->size + 1, 1.3);
            for (i = 0; i < cuts->count; i++)
                cuts->cuts[i]->cliques = cuts->cliques->cliques;
        }
        y                         = cuts->cliques->size++;
        cuts->cliques->cliques[y] = clique = malloc(sizeof(graph_clique));
    }

    clique_copy(c, clique);

    clique->refcount             = 1;
    clique->hashnext             = cuts->cliques->cliquehash[x];
    cuts->cliques->cliquehash[x] = y;

    return y;
}

static void
cp_unregister_cut_repo_clique(cp_cut_repo *cuts, int c)
{
    int x, y, yprev;

    cuts->cliques->cliques[c]->refcount--;
    if (cuts->cliques->cliques[c]->refcount)
        return;
    x = clique_hash(cuts->cliques->cliques[c]) % cuts->cliques->cliquehashsize;
    y = cuts->cliques->cliquehash[x];
    if (y == c)
    {
        cuts->cliques->cliquehash[x] = cuts->cliques->cliques[c]->hashnext;
    }
    else
    {
        yprev = y;
        y     = cuts->cliques->cliques[y]->hashnext;
        while (y != c && y != -1)
        {
            yprev = y;
            y     = cuts->cliques->cliques[y]->hashnext;
        }
        if (y == -1)
        {
            fprintf(stderr, "Couldn't find clique to delete from hash\n");
            return;
        }
        cuts->cliques->cliques[yprev]->hashnext =
        cuts->cliques->cliques[c]->hashnext;
    }
    free(cuts->cliques->cliques[c]->nodes);
    cuts->cliques->cliques[c]->nodes    = NULL;
    cuts->cliques->cliques[c]->segcount = -1;
    cuts->cliques->cliques[c]->hashnext = cuts->cliques->cliquefree;
    cuts->cliques->cliquefree           = c;
}

static int
sort_cliques(const void *xx, const void *yy)
/**************************************************************************/
{
    int x = *(int *)xx;
    int y = *(int *)yy;

    if (x < y)
        return -1;
    if (x > y)
        return +1;

    return 0;
}

int
cp_register_cut_repo_cliques(cp_cut_repo *cuts, cp_cut *cut)
{
    int i, j;

    if (cut->tcount)
    {
        cut->teeth_cid = malloc(cut->tcount * sizeof(int));
        for (i = 0; i < cut->tcount; i++)
        {
            cut->teeth_cid[i] =
            cp_register_cut_repo_clique(cuts, cut->teeth[i]);
            clique_free(&cut->teeth[i]);
            if (cut->teeth_cid[i] == -1)
            {
                for (j = 0; j < i; j++)
                    cp_unregister_cut_repo_clique(cuts, cut->teeth_cid[j]);
                free(cut->teeth_cid);
                cut->cliques = NULL;
                return 1;
            }
        }
        free(cut->teeth);
        cut->teeth = NULL;
    }
    if (cut->hcount)
    {
        cut->handle_cid = malloc(cut->hcount * sizeof(int));
        for (i = 0; i < cut->hcount; i++)
        {
            cut->handle_cid[i] =
            cp_register_cut_repo_clique(cuts, cut->handles[i]);
            clique_free(&cut->handles[i]);
            free(cut->handles[i]);

            if (cut->handle_cid[i] == -1)
            {
                for (j = 0; j < i; j++)
                    cp_unregister_cut_repo_clique(cuts, cut->handle_cid[j]);
                free(cut->handle_cid);
                cut->cliques = NULL;
                return 1;
            }
        }
        free(cut->handles);
        cut->handles = NULL;
    }
    cut->cliques = cuts->cliques->cliques;
    cut->age     = 0;

    qsort(cut->handle_cid, cut->hcount, sizeof(int), sort_cliques);
    qsort(cut->teeth_cid, cut->tcount, sizeof(int), sort_cliques);

    return 0;
}

int
cp_unregister_cut_repo_cliques(cp_cut_repo *cuts, cp_cut *cut)
{
    int i;

    if (cut->hcount + cut->tcount)
    {
        assert(cuts->cliques->cliques == cut->cliques);
        if (cut->hcount)
        {
            cut->handles = malloc(cut->hcount * sizeof(graph_clique *));
            for (i = 0; i < cut->hcount; i++)
            {
                cut->handles[i] = clique_create();
                clique_copy(cuts->cliques->cliques[cut->handle_cid[i]],
                            cut->handles[i]);
                cp_unregister_cut_repo_clique(cuts, cut->handle_cid[i]);
            }
            free(cut->handle_cid);
            cut->handle_cid = NULL;
        }
        if (cut->tcount)
        {
            cut->teeth = malloc(cut->tcount * sizeof(graph_clique *));
            for (i = 0; i < cut->tcount; i++)
            {
                cut->teeth[i] = clique_create();
                clique_copy(cuts->cliques->cliques[cut->teeth_cid[i]],
                            cut->teeth[i]);
                cp_unregister_cut_repo_clique(cuts, cut->teeth_cid[i]);
            }
            free(cut->teeth_cid);
            cut->teeth_cid = NULL;
        }
        cut->cliques = NULL;
        cut->age     = 0;
    }

    return 0;
}
