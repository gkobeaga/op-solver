#ifndef CP_EXACT_BAC_H
#define CP_EXACT_BAC_H

#include "../../../ip/ip.h"
#include "../exact.h"

typedef struct cp_exact_bac_param
{
    long time_limit;
    int sep_logical;
    int sep_sec_comps;
    int sep_sec_exact;
#define CP_SEC_NONE (0)
#define CP_SEC_HONG_ST (1)
#define CP_SEC_HONG_TS (2)
#define CP_SEC_HONG (CP_SEC_HONG_ST)
#define CP_SEC_HONG_BOTH (CP_SEC_HONG_ST | CP_SEC_HONG_TS)
#define CP_SEC_GOMORYHU (4)
    int sep_sec_cc_2;
    int sep_sec_cc_extra;
    int sep_clique_heur;
    int sep_clique_exact;
    int sep_clique_mincut;
    int sep_blossom_mst;
    int sep_blossom_fast;
    int sep_blossom_ghfast;
    int sep_cover_edge;
    int sep_cover_vertex;
    int sep_cover_cycle;
    int sep_path;
    int sep_loop;
#define CP_SEP_LOOP_BI_LEVEL (0)
#define CP_SEP_LOOP_THREE_LEVEL (1)
    int srk_rule;
#define CP_SRK_NONE (0)
#define CP_SRK_C1 (1)
#define CP_SRK_C1C2 (2)
#define CP_SRK_C1C2C3 (3)
#define CP_SRK_S1 (4)
    int srk_s2;
    int srk_s3;
    int srk_extra;
    int sec_max_vout;
    int sec_max_vin;
    int sec_max_viol;
    int sec_max_cut;
    int sec_max_cut_x_clique;
    int xheur_vph;
    int xheur_vph_meta;
    double pruning_tol;
} cp_exact_bac_param;

typedef struct cp_exact_bac_stats
{
    stats_item *total;
    stats_item *sep_loop;
    stats_item *sep_loop_it;
    stats_item *sep_loop_inner;
    stats_item *sep_loop_inner_it;
    stats_item *sep_loop_middle;
    stats_item *sep_loop_middle_it;
    stats_item *sep_loop_outer;
    stats_item *sep_loop_outer_it;
    stats_item *sep_logical;
    stats_item *sep_sec_comps;
    stats_item *sep_sec_exact;
    stats_item *sep_blossom_fast;
    stats_item *sep_blossom_ghfast;
    stats_item *sep_blossom_mst;
    stats_item *sep_connect_mincut;
    stats_item *sep_cover_edge;
    stats_item *sep_cover_vertex;
    stats_item *sep_cover_cycle;
    stats_item *sep_path;
    stats_item *age_cuts;
    stats_item *age_vars;
    stats_item *add_cuts;
    stats_item *add_vars;
    stats_item *xheur_branch;
    stats_item *xheur_sep;
    stats_item *lp_opt;
    stats_item *misc;
    const char *file;
} cp_exact_bac_stats;

typedef struct kp_env kp_env;
typedef struct cp_cut cp_cut;
typedef struct cp_cut_repo cp_cut_repo;
typedef struct ip_exact_bac_env ip_exact_bac_env;
typedef struct cp_exact_bac_env
{
    int verbosity;
    cp_exact_bac_param *param;
    cp_exact_bac_stats *stats;
    cp_init_env *init;
    cp_heur_env *heur;
    ip_exact_bac_env *ip;
    cp_cut_repo *cuts;
    cp_cut *cut_queue;
    kp_env *kp;
    int ngot;
} cp_exact_bac_env;

cp_exact_bac_stats *
cp_create_exact_bac_stats(void);
void
cp_free_exact_bac_stats(cp_exact_bac_stats **stats);
int
cp_write_exact_bac_stats(cp_prob *cp, cp_exact_bac_env *env);

cp_exact_bac_param *
cp_create_exact_bac_param(void);
void
cp_free_exact_bac_param(cp_exact_bac_param **param);
int
cp_parse_exact_bac_args(int argc, char *argv[], cp_exact_bac_env *env);

cp_exact_bac_env *
cp_create_exact_bac_env(void);
void
cp_free_exact_bac_env(cp_exact_bac_env **env);

typedef struct cp_cut cp_cut;
typedef struct cp_cut_logic cp_cut_logic;
typedef struct cp_cut_cover_edge cp_cut_cover_edge;
typedef struct cp_cut_cover_vertex cp_cut_cover_vertex;
typedef struct cp_cut_connect cp_cut_connect;
typedef struct cp_cut_path cp_cut_path;

typedef struct cp_cut_skeleton
{
    int atomcount;
    int *atoms;
    int nverts;
    int *verts;
} cp_cut_skeleton;

typedef struct cp_cut
{
    int hcount;
    int tcount;

    int *handle_cid;
    graph_clique **handles;
    int *teeth_cid;
    graph_clique **teeth;

    int *verts;
    double vycoef;
    graph_clique **cliques; // Pointer to a repo cliques
    cp_cut_logic *logical;
    cp_cut_cover_edge *cover_edge;
    cp_cut_cover_vertex *cover_vertex;
    cp_cut_connect *connect;
    cp_cut_path *path;
    int min_depth;

    int age;
    double rhs;
    char sense;
    char branch;

    cp_cut *next;
    cp_cut *prev;
    cp_cut *hash_prev;
    cp_cut *hash_next;

    cp_cut_skeleton *skel;

} cp_cut;

cp_cut *
cp_create_cut(void);
void
cp_free_cut(cp_cut **cut);
double
cp_eval_cut(solver_graph *graph, cp_cut *cut);
void
cp_get_cut_arcs(solver_graph *graph, cp_cut *cut, int *nzcnt,
                graph_arc **nzlist);
void
cp_print_cut(cp_cut *cut);
int
cp_clear_depth_invalid_cuts(cp_prob *cp, cp_exact_bac_env *bac_env, int depth,
                            int *ndel);

typedef struct cp_cut_logic
{
    int arc[2];
    graph_arc *arcptr;
    int v;
} cp_cut_logic;

typedef struct cp_cut_cover_edge
{
    int *arcs;
    int na;
    int strong; // Strong if cycle cover: x(T) <= Y(V(T)) - 1
} cp_cut_cover_edge;

typedef struct cp_cut_cover_vertex
{
    graph_clique *verts;
    double rhs;
    char sense;
} cp_cut_cover_vertex;

typedef struct cp_cut_connect
{
    graph_clique *verts;
    double rhs;
    char sense;
} cp_cut_connect;

typedef struct cp_cut_path
{
    int *arcs;
    int na;
    int *farcs;
    int fna;
} cp_cut_path;

typedef struct cp_cut_hash cp_cut_hash;
typedef struct cp_cut_repo
{
    int count;
    cp_cut **cuts;
    int space;
    cp_cut_hash *hash;
    graph_clique_repo *cliques;
} cp_cut_repo;

typedef struct cp_cut_hash
{
    int nelem;
    int maxelem;
    int size;
    double maxdensity;
    double lowdensity;
    cp_cut **table;
} cp_cut_hash;

cp_cut_hash *
cp_create_cut_hash(int size);
unsigned int
cp_get_cut_hash(cp_cut *cut);
int
cp_eq_cut(cp_cut *cut1, cp_cut *cut2),
cp_add_cut_hash(cp_cut_hash *hash, cp_cut *cut),
cp_del_cut_hash(cp_cut_hash *hash, cp_cut *cut);
cp_cut *
cp_find_cut_hash(cp_cut_hash *hash, cp_cut *cut);

cp_cut_skeleton *
cp_create_cut_skeleton(void);

typedef struct cp_cut_list_item
{
    cp_cut *cut;
    double val;
} cp_cut_list_item;

typedef struct cp_cut_list
{
    cp_cut_list_item *cuts;
    int max_size;
    int size;
    double max_val;
} cp_cut_list;

cp_cut_list *
cp_create_cut_list(int count, double max_val);
void
cp_insert_cut_list(cp_cut_list *list, double val, cp_cut *cut, int *finish),
cp_append_cut_list(cp_cut_list *main_list, cp_cut_list *from_list, int *stop),
cp_get_cut_list(cp_prob *cp, cp_cut_list *list, int *cutcount, cp_cut **cuts),
cp_erase_cut_list(cp_cut_list *cut_list), cp_free_cut_list(cp_cut_list **list);

int
cp_opt_exact_bac(cp_prob *op, cp_exact_bac_env *env, cp_sol *sol);
int
cp_init_exact_bac(cp_prob *op, cp_exact_bac_env *env);

void
cp_free_cut_skeleton(cp_cut_skeleton **skel);
int
cp_copy_cut_skeleton(cp_cut_skeleton *in, cp_cut_skeleton *out);
int
cp_build_cut_skeleton(cp_cut *cut, int nodecount);
int
cp_eq_cut_skeletons(cp_cut_skeleton *a, cp_cut_skeleton *b);

cp_cut_repo *
cp_create_cut_repo(cp_prob *cp);
void
cp_free_cut_repo(cp_cut_repo **repo);
int
cp_register_cut_repo_cliques(cp_cut_repo *cuts, cp_cut *cut),
cp_unregister_cut_repo_cliques(cp_cut_repo *cuts, cp_cut *cut);

int
cp_sep_loop(void *prob, void *env, int *nadded_tot);

int
cp_sep_logical(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
               cp_cut **cuts),
cp_sep_sec_exact(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                 cp_cut **cuts),
cp_sep_sec_exact_hong(cp_prob *cp, cp_exact_bac_env *bac_env,
                      solver_graph *graph, graph_clique_repo *repo),
cp_sep_sec_exact_gomoryhu(cp_prob *cp, cp_exact_bac_env *bac_env,
                          solver_graph *graph, graph_clique_repo *repo),
cp_sep_sec_comps(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                 cp_cut **cuts),
cp_find_sec_cliques(solver_graph *srkgraph, graph_vertex **sverts, int scount,
                    double cutval, double lowerbound, graph_clique_repo *repo),
cp_sep_blossom_fast(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                    cp_cut **cuts),
cp_sep_blossom_ghfast(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                      cp_cut **cuts),
cp_sep_blossom_mst(cp_prob *cp, cp_exact_bac_env *env, int *cutcount,
                   cp_cut **cuts),
cp_sep_cover_edge(cp_prob *cp, cp_exact_bac_env *env, int *cutcount,
                  cp_cut **cuts),
cp_sep_cover_vertex(cp_prob *cp, cp_exact_bac_env *env, int *cutcount,
                    cp_cut **cuts),
cp_sep_cover_cycle(cp_prob *cp, cp_exact_bac_env *env, int *cutcount,
                   cp_cut **cuts),
cp_sep_path(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
            cp_cut **cuts);

int
cp_get_clique_repo_sec_cuts(cp_prob *cp, cp_exact_bac_env *bac_env,
                            solver_graph *graph, graph_clique_repo *repo,
                            cp_cut_list *cut_list);

cp_cut *
cp_conv_cut_sol2connect(cp_sol *sol),
*cp_conv_cut_verts2connect(solver_graph *graph, int vcount,
                           graph_vertex **verts);

void
cp_get_cut_arcs(solver_graph *graph, cp_cut *cut, int *nzcnt,
                graph_arc **nzlist);

int
cp_shrink_exact_bac_graph(cp_prob *cp, cp_exact_bac_env *bac_env,
                          solver_graph *graph, graph_vertex *qstart,
                          graph_clique_repo *repo),
cp_shrink_exact_bac_graph_c1(cp_prob *cp, cp_exact_bac_env *bac_env,
                             solver_graph *graph, graph_vertex *qstart,
                             graph_clique_repo *repo),
cp_shrink_exact_bac_graph_c1c2(cp_prob *cp, cp_exact_bac_env *bac_env,
                               solver_graph *graph, graph_vertex *qstart,
                               graph_clique_repo *repo),
cp_shrink_exact_bac_graph_c1c2c3(cp_prob *cp, cp_exact_bac_env *bac_env,
                                 solver_graph *graph, graph_vertex *qstart,
                                 graph_clique_repo *repo),
cp_shrink_exact_bac_graph_s1(cp_prob *cp, cp_exact_bac_env *bac_env,
                             solver_graph *graph, graph_vertex *qstart,
                             graph_clique_repo *repo),
cp_shrink_exact_bac_graph_s1s2(cp_prob *cp, cp_exact_bac_env *bac_env,
                               solver_graph *graph, graph_vertex *qstart,
                               graph_clique_repo *repo);

#define ADD_TO_SRK_QUEUE(n)                                                    \
    {                                                                          \
        if (!(n)->onqueue)                                                     \
        {                                                                      \
            (n)->qnext = NULL;                                                 \
            if (qtail)                                                         \
                qtail->qnext = (n);                                            \
            else                                                               \
                qhead = (n);                                                   \
            qtail        = (n);                                                \
            (n)->onqueue = 1;                                                  \
        }                                                                      \
    }

int
cp_recover_infeas(void *prob, void *env),
cp_add_badvars(void *prob, void *env, int *nadded);

typedef struct ip_branch ip_branch;
int
cp_find_branch_edge(cp_prob *cp, cp_exact_bac_env *bac_env, graph_arc **edge),
cp_find_branch(void *prob, void *env, ip_branch **bobj);

int
cp_get_branch_dual_bound(void *prob, void *env, mpf_t bound),
cp_verify_branch_infeas(void *prob, void *env, int *yesno),
cp_verify_branch_prune(void *prob, void *env, int *yesno);

int
cp_age_cuts(cp_prob *cp, cp_exact_bac_env *bac_env, int *ndel),
cp_age_vars(cp_prob *cp, cp_exact_bac_env *bac_env, int *ndel);

/* Lib: LP */
typedef struct lp_data lp_data;
lp_data *
cp_build_lp_data(cp_prob *cp, cp_exact_bac_env *bac_env, int col_start,
                 int col_end);
int
cp_update_lp(cp_prob *cp, cp_exact_bac_env *bac_env, int age),
cp_update_lp_sol(cp_prob *cp, cp_exact_bac_env *bac_env);

int
cp_add_lp_arcs(cp_prob *cp, cp_exact_bac_env *bac_env, graph_arc **prlist,
               int narcs),
cp_add_lp_cuts(cp_prob *cp, cp_exact_bac_env *bac_env, cp_cut **cuts,
               int *nadded),
cp_add_lp_cut2data(cp_prob *cp, cp_exact_bac_env *bac_env, cp_cut *cut,
                   lp_data **data);

int
cp_add_lp_cut(cp_prob *prob, cp_exact_bac_env *bac_env, cp_cut *cut),
cp_del_lp_cut(cp_prob *prob, cp_exact_bac_env *bac_env, int ind);

int
cp_get_xheur_greedy__(void *cp_, void *bac_env_, void *sol_, int *nadded),
cp_get_xheur_greedy(cp_prob *cp, cp_exact_bac_env *bac_env, cp_sol *sol),
cp_get_xheur_ea__(void *cp_, void *bac_env_, void *sol_, int *nadded),
cp_get_xheur_ea(cp_prob *cp, cp_exact_bac_env *bac_env, cp_sol *sol);

#define inside_pos(n, i, j)                                                    \
    (i < j ? (i * (((n)-1) + ((n)-i)) / 2 + j - i - 1)                         \
           : (j * (((n)-1) + ((n)-j)) / 2 + i - j - 1))

int
cp_is_sol_lp_feasible(cp_prob *cp, cp_exact_bac_env *bac_env, cp_sol *sol);
#endif
