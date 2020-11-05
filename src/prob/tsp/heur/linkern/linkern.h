#ifndef TSP_HEUR_LINKERN_H
#define TSP_HEUR_LINKERN_H

typedef struct tsp_heur_linkern_param
{
    long time_limit;
    double length_bound;
    int nruns; /* Number of repetitions */
    int kick_type;
#define TSP_LINKERN_RANDOM_KICK (0)
#define TSP_LINKERN_GEOMETRIC_KICK (1)
#define TSP_LINKERN_CLOSE_KICK (2)
#define TSP_LINKERN_WALK_KICK (3)
    int nkicks;
} tsp_heur_linkern_param;

typedef struct stats_item stats_item;
typedef struct tsp_heur_linkern_stats
{
    stats_item *total;
    stats_item *steps;
} tsp_heur_linkern_stats;

typedef struct tsp_heur_linkern_env
{
    int verbosity;
    tsp_heur_linkern_param *param;
    tsp_heur_linkern_stats *stats;
    int ecount;
    int *elist;
} tsp_heur_linkern_env;

tsp_heur_linkern_stats *
tsp_create_heur_linkern_stats(void);
void
tsp_free_heur_linkern_stats(tsp_heur_linkern_stats **stats);

tsp_heur_linkern_param *
tsp_create_heur_linkern_param(void);
void
tsp_free_heur_linkern_param(tsp_heur_linkern_param **param);

tsp_heur_linkern_env *
tsp_create_heur_linkern_env(void);
void
tsp_free_heur_linkern_env(tsp_heur_linkern_env **env);

typedef struct cp_prob tsp_prob;
typedef struct cp_sol tsp_sol;
int
tsp_opt_heur_linkern(tsp_prob *tsp, tsp_heur_linkern_env *env, tsp_sol *sol);

/****************************************************************************/
/*                                                                          */
/*                  FLIPPER HEADER (TWO-LIST)                               */
/*                                                                          */
/****************************************************************************/

typedef struct __lk_parentnode __lk_parentnode;
typedef struct __lk_childnode __lk_childnode;
typedef struct __lk_parentnode
{
    __lk_parentnode *adj[2];
    __lk_childnode *ends[2];
    int size;
    int id;
    int rev;
} __lk_parentnode;

typedef struct __lk_childnode
{
    __lk_parentnode *parent;
    __lk_childnode *adj[2];
    int id;
    int name;
} __lk_childnode;

typedef struct __lk_flipper
{
    __lk_parentnode *parents;
    __lk_childnode *children;
    int reversed;
    int nsegments;
    int groupsize;
    int split_cutoff;
} __lk_flipper;

typedef struct solver_data solver_data;
int
__linkern_flipper_init(__lk_flipper *f, tsp_sol *sol),
__linkern_flipper_next(__lk_flipper *f, int x),
__linkern_flipper_prev(__lk_flipper *f, int x),
__linkern_flipper_sequence(__lk_flipper *f, int x, int y, int z);
void
__linkern_flipper_flip(__lk_flipper *F, int x, int y),
__linkern_flipper_cycle(__lk_flipper *F, tsp_sol *sol, solver_data *data),
__linkern_flipper_finish(__lk_flipper *F);

#endif /* __FLIPPER_H */
