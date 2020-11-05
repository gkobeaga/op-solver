#include "op-solver.h"

graph_clique_repo *
clique_create_repo(int n)
{
    int i;
    graph_clique_repo *repo = malloc(sizeof(graph_clique_repo));
    repo->size              = 0;
    repo->cliquespace       = 1000;
    repo->cliquehashsize    = 0;
    repo->cliquehash        = NULL;
    repo->cliques = malloc(repo->cliquespace * sizeof(graph_clique *));
    repo->tot_n   = n;

    repo->cliquehashsize = prime_next((unsigned int)2 * n);
    repo->cliquehash     = malloc(repo->cliquehashsize * sizeof(int));
    if (!repo->cliquehash)
    {
        repo->cliquehashsize = 0;
        return NULL;
    }
    for (i = 0; i < repo->cliquehashsize; i++)
    {
        repo->cliquehash[i] = -1;
    }
    repo->cliquefree = -1;
    return repo;
}

void
clique_free_repo(graph_clique_repo **repo)
{
    int i;
    if (*repo)
    {
        if ((*repo)->cliquehash)
            free((*repo)->cliquehash);
        if ((*repo)->cliques)
        {
            for (i = 0; i < (*repo)->size; i++)
            {
                clique_free(&(*repo)->cliques[i]);
            }
            free((*repo)->cliques);
        }
        free(*repo);
    }
}

int
clique_register_repo(solver_graph *graph, graph_clique_repo *repo, graph_clique *c)
{
    int x = clique_hash(c) % repo->cliquehashsize;
    int y = repo->cliquehash[x];
    graph_clique *clique;

    while (y != -1)
    {
        clique = repo->cliques[y];
        if (clique_eq(c, clique))
        {
            clique->refcount++;
            return y;
        }
        y = clique->hashnext;
    }

    if (repo->cliquefree != -1)
    {
        y                = repo->cliquefree;
        clique           = repo->cliques[y];
        repo->cliquefree = clique->hashnext;
        clique           = repo->cliques[y];
    }
    else
    {
        if (repo->size + 1 >= repo->cliquespace)
        {
            realloc_scale(repo->cliques, repo->cliquespace, repo->size + 1,
                          1.3);
        }
        y                = repo->size++;
        repo->cliques[y] = clique = malloc(sizeof(graph_clique));
    }

    clique_copy(c, clique);

    clique->refcount    = 1;
    clique->hashnext    = repo->cliquehash[x];
    repo->cliquehash[x] = y;

    return y;
}

void
clique_unregister_repo(solver_graph *graph, graph_clique_repo *repo, int c)
{
    int x, y, yprev;

    repo->cliques[c]->refcount--;
    if (repo->cliques[c]->refcount)
        return;
    x = clique_hash(repo->cliques[c]) % repo->cliquehashsize;
    y = repo->cliquehash[x];
    if (y == c)
    {
        repo->cliquehash[x] = repo->cliques[c]->hashnext;
    }
    else
    {
        yprev = y;
        y     = repo->cliques[y]->hashnext;
        while (y != c && y != -1)
        {
            yprev = y;
            y     = repo->cliques[y]->hashnext;
        }
        if (y == -1)
        {
            fprintf(stderr, "Couldn't find clique to delete from hash\n");
            return;
        }
        repo->cliques[yprev]->hashnext = repo->cliques[c]->hashnext;
    }
    free(repo->cliques[c]->nodes);
    repo->cliques[c]->segcount = -1;
    repo->cliques[c]->hashnext = repo->cliquefree;
    repo->cliquefree           = c;
}

int
clique_register_repo_srkvertices(solver_graph *graph, graph_vertex *u,
                                 graph_vertex *v, double eweight,
                                 graph_clique_repo *repo)
{
    int rval = 0;

    int i, t1vcount, t2vcount, cvcount;
    graph_vertex **t1verts = NULL, **t2verts = NULL, **cverts = NULL;
    graph_clique *clique = NULL;

    rval = graph_expand_vertex(graph, u, &t1vcount, &t1verts);
    check_rval(rval, "graph_expand_vertex failed", CLEANUP);
    rval = graph_expand_vertex(graph, v, &t2vcount, &t2verts);
    check_rval(rval, "graph_expand_vertex failed", CLEANUP);

    cvcount = t1vcount + t2vcount;
    cverts  = malloc(cvcount * sizeof(graph_vertex *));

    for (i = 0; i < t1vcount; i++) cverts[i] = t1verts[i];
    for (i = 0; i < t2vcount; i++) cverts[t1vcount + i] = t2verts[i];

    if (cvcount <= graph->nv / 2 && cvcount > 2 && cvcount < graph->nv - 2)
        clique = clique_conv_vertices2clique(graph->orig, cverts, cvcount);
    else
        clique = clique_conv_vertices2coclique(graph->orig, cverts, cvcount);
    clique->val = 2 * u->y + 2 * v->y - 2 * eweight;
    clique_register_repo(graph->orig, repo, clique);

CLEANUP:
    if (clique)
        clique_free(&clique);
    if (t1verts)
        free(t1verts);
    if (t2verts)
        free(t2verts);
    if (cverts)
        free(cverts);
    return rval;
}
