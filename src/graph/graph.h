#ifndef GRAPH_H
#define GRAPH_H

#define NV_MAX 1000000
#define NA_MAX 6000000

typedef struct solver_graph solver_graph;
typedef struct graph_vertex graph_vertex;
typedef struct graph_arc graph_arc;

typedef struct graph_vertex
{
    int i;
    int ind;
    int fixed;
    int branch;
    int deg;
    int mark;
    double obj;
    double coef;
    double y;
    double rc;
    int age;
    int used;
    graph_vertex *shrunk;

    int comp;
    int depth;
    int lowdepth;
    graph_arc *in;
    graph_arc *out;
    graph_arc *edge;

    graph_vertex *next;
    graph_vertex *prev;

    graph_vertex *fixed_next;
    graph_vertex *fixed_prev;

    // SRK
    graph_vertex *members;
    int nmembers;
    graph_vertex *orig;
    graph_vertex *max;

    // SRK rules (C2, C3)
    graph_vertex *qnext;
    double vt_x;
    double ut_x;

    // Maxflow
    double excess;
    int active;
    int label;
    graph_vertex *search_next;
    graph_vertex *high_next;
    graph_vertex *level_prev;
    graph_vertex *level_next;
    graph_arc *ecurrent;

    int onecnt; // TODO
    int onqueue;
    int mark_aux;

    graph_vertex *parent;
} graph_vertex;

typedef struct graph_arc
{
    graph_vertex *tail;
    graph_vertex *head;

    graph_arc *prev;
    graph_arc *next;

    graph_arc *t_prev;
    graph_arc *t_next;
    graph_arc *h_prev;
    graph_arc *h_next;
    graph_arc *rev;
    graph_arc *orig;

    graph_arc *hash_prev;
    graph_arc *hash_next;

    int i;
    int ind;
    int age;
    int aux;

    double obj;
    double cost;

    int fixed;
    int branch;

    double x;
    double rc;
    double coef;
    int used;

    double flow;
} graph_arc;

typedef struct graph_arc_hash
{
    graph_arc **table;
    unsigned int size;
    unsigned int mult;
} graph_arc_hash;

typedef struct solver_data solver_data;
typedef struct solver_graph
{
    solver_data *data;
    int nv_space;
    int nv;
    int na_space;
    int na;
    int directed;
    int marker;
    int connected;
    // int integral;
    solver_graph *shrunk;
    graph_arc_hash *archash;
    graph_vertex **v;
    graph_vertex *fixed;
    graph_arc **arcs;
    solver_graph *orig;

    // Number of non-null vertices
    int n3v;

    graph_vertex *tail;
    graph_vertex *head;

    int onecnt;
    int original_ncount;
    int original_ecount;
} solver_graph;

solver_graph *
graph_create(void);
void
graph_erase(solver_graph *graph),
graph_free(solver_graph **graph);
int
graph_copy(solver_graph *ingraph, solver_graph *outgraph);

/* Used to break ties when reordering
 * 0 Random (Much slower than the others, around 3 times,
 *           specially if the shrinking is not considered)
 * 1 Maintain graph order
 * 2 Maintain previous order  */
#define REORDER 0

int
graph_add_vertices(solver_graph *graph, int nadd),
graph_reorder_vertices(solver_graph *graph);
void
graph_del_vertices(solver_graph *graph, int ndel, const int num[]);

graph_arc *
graph_add_arc(solver_graph *graph, int i, int j);
graph_arc *
graph_find_arc(solver_graph *graph, graph_vertex *tail, graph_vertex *head);
void
graph_del_arc(solver_graph *graph, graph_arc **a);

void
graph_identify_vertices(solver_graph *graph, graph_vertex *v, graph_vertex *u);
int
graph_expand_vertex(solver_graph *graph, graph_vertex *srkv, int *vcount,
                    graph_vertex ***verts);

int
graph_init_arc_hash(solver_graph *graph, int size),
graph_del_arc_hash(graph_arc_hash *hash, graph_arc *arc),
graph_getall_arc_hash(graph_arc_hash *hash, int *narcs, graph_arc ***arcs);
void
graph_add_arc_hash(graph_arc_hash *hash, graph_arc *arc),
graph_free_arc_hash(graph_arc_hash **hash);
graph_arc *
graph_find_arc_hash(graph_arc_hash *hash, int end0, int end1);

#define otherend(arc, v) ((arc)->tail->i == (v->i) ? (arc)->head : (arc)->tail)
#define outnext(arc, v)                                                        \
    ((arc)->tail->i == (v->i) ? (arc)->t_next : (arc)->h_next)
#define outprev(arc, v)                                                        \
    ((arc)->tail->i == (v->i) ? (arc)->t_prev : (arc)->h_prev)
#define innext(arc, v)                                                         \
    ((arc)->head->i == (v->i) ? (arc)->t_next : (arc)->h_next)
#define inprev(arc, v)                                                         \
    ((arc)->head->i == (v->i) ? (arc)->t_prev : (arc)->h_prev)

#define set_tnext(e, v, next)                                                  \
    if ((e)->tail->i == (v)->i)                                                \
    {                                                                          \
        (e)->t_next = next;                                                    \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        (e)->h_next = next;                                                    \
    }

#define set_tprev(e, v, prev)                                                  \
    if ((e)->tail->i == (v)->i)                                                \
    {                                                                          \
        (e)->t_prev = prev;                                                    \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        (e)->h_prev = prev;                                                    \
    }

#define set_hnext(e, v, next)                                                  \
    if ((e)->head->i == (v)->i)                                                \
    {                                                                          \
        (e)->h_next = next;                                                    \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        (e)->t_next = next;                                                    \
    }

#define set_hprev(e, v, prev)                                                  \
    if ((e)->head->i == (v)->i)                                                \
    {                                                                          \
        (e)->h_prev = prev;                                                    \
    }                                                                          \
    else                                                                       \
    {                                                                          \
        (e)->t_prev = prev;                                                    \
    }

typedef struct graph_clique_segment
{
    int lo;
    int hi;
} graph_clique_segment;

typedef struct graph_clique graph_clique;
typedef struct graph_clique
{
    int i;
    int nodecount;
    int segcount;
    graph_clique_segment *nodes;
    int hashnext;
    int refcount;
    double val;
    graph_clique *prev;
    graph_clique *next;
    int tot_n;
} graph_clique;

#define FOREACH_NODE_IN_CLIQUE(i, c, tmp)                                      \
    for (tmp = 0; tmp < c->segcount; tmp++)                                    \
        for (i = c->nodes[tmp].lo; i <= c->nodes[tmp].hi; i++)

static inline void
__foreach_not_in_next_tmp(graph_clique *c, int *nseg, int *i)
{
    (*i)++;
    if (*nseg < c->segcount)
    {
        int tmp = *nseg;
        *nseg   = *i >= c->nodes[*nseg].lo ? (*nseg) + 1 : *nseg;
        for (; *i >= c->nodes[tmp].lo && *i <= c->nodes[tmp].hi; (*i)++);
    }
}

#define FOREACH_NODE_NOT_IN_CLIQUE(c, i)                                       \
    static int nseg; /* TODO:Is this safe in parallel? */                      \
    for (nseg = 0, i = -1, __foreach_not_in_next_tmp((c), &(nseg), &(i));      \
         i < (c)->tot_n; __foreach_not_in_next_tmp((c), &(nseg), &(i)))

int
graph_get_comps(solver_graph *graph, int *ncomp, int **compscount,
                graph_vertex ***comps);

graph_clique *
clique_create(void),
*clique_conv_vertices2clique(solver_graph *graph, graph_vertex **vert,
                             int vcount),
*clique_conv_vertices2coclique(solver_graph *graph, graph_vertex **vert,
                               int vcount);

int
clique_conv_array2clique(int *vert, int vcount, graph_clique **clique),
clique_conv_clique2array(graph_clique *clique, int **ar, int *count);
int
clique_copy(graph_clique *in, graph_clique *out),
clique_eq(graph_clique *c, graph_clique *d), clique_count(graph_clique *clique);
void
clique_free(graph_clique **clique),
clique_print(graph_clique *clique);
unsigned int
clique_hash(graph_clique *c);

typedef struct graph_clique_repo
{
    int size;
    int cliquespace;
    int cliquehashsize;
    int cliquefree;
    int *cliquehash;
    graph_clique **cliques;
    int tot_n;
} graph_clique_repo;

graph_clique_repo *
clique_create_repo(int n);
int
clique_register_repo(solver_graph *graph, graph_clique_repo *repo,
                     graph_clique *c);
void
clique_unregister_repo(solver_graph *graph, graph_clique_repo *repo, int c),
clique_free_repo(graph_clique_repo **repo);

int
clique_register_repo_srkvertices(solver_graph *graph, graph_vertex *u,
                                 graph_vertex *v, double eweight,
                                 graph_clique_repo *repo);

void
graph_print(solver_graph *graph),
graph_write(solver_graph *graph, char *fname);
solver_graph *
graph_read(char *fname);

/* Minimum and Maximum Spanning trees*/

#define SOLVER_MST_MAX 1
#define SOLVER_MST_MIN 0
solver_graph *
graph_get_mst(solver_graph *graph, int max),
*graph_get_mst_max(solver_graph *graph),
*graph_get_mst_min(solver_graph *graph);

/* Minimum and Maximum flow problems*/
double
graph_get_maxflow_st(solver_graph *graph, graph_vertex *s, graph_vertex *t);
int
graph_get_mincut_st(solver_graph *graph, graph_vertex *s, graph_vertex *t,
                    double *value, graph_vertex ***verts, int *vcount);

typedef struct graph_ghtree_node graph_ghtree_node;
typedef struct graph_ghtree_node
{
    int i;
    int mark;
    graph_ghtree_node *parent;
    graph_ghtree_node *next_sibling;
    graph_ghtree_node *prev_sibling;
    graph_ghtree_node *child;
    double cutval;
    int nchild;
    int ndesc;
    graph_vertex *in_max;
    graph_vertex *out_max;
    graph_vertex *special;
    int mcount;
    graph_ghtree_node **members;
    graph_vertex *pseudonode;
    graph_ghtree_node *next;
} graph_ghtree_node;

typedef struct graph_ghtree
{
    int nn;
    graph_ghtree_node **n;
    graph_ghtree_node *root;
    int marker;
} graph_ghtree;

graph_ghtree *
graph_create_ghtree(void);
graph_ghtree_node *
graph_add_ghtree_node(graph_ghtree *ghtree);
void
graph_add_ghtree_nodes(graph_ghtree *ghtree, int nn),
graph_print_ghtree(graph_ghtree *ghtree),
graph_free_ghtree(graph_ghtree *ghtree);

int
graph_get_ghtree(solver_graph *graph, int rootind, graph_ghtree **ghtree);

void
graph_plot(solver_graph *graph);
#endif
