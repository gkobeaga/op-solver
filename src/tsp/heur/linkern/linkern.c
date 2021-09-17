#include "op-solver.h"
#include "tsp/tsp.h"
#include "tsp/heur/heur.h"
#include "linkern.h"

#define MAXDEPTH 25 /* Shouldn't really be less than 2.             */
#define KICK_MAXDEPTH 50

#define MARK_LEVEL 10 /* Number of tour neighbors after 4-swap kick   */
#define BACKTRACK 4
#define MAX_BACK 12 /* Upper bound on the XXX_count entries         */
static const int backtrack_count[BACKTRACK] = {4, 3, 3, 2};
static const int weird_backtrack_count[3]   = {4, 3, 3};

#define Edgelen(n1, n2, data) data_get_norm(data, n1, n2)

#define FLIP(aprev, a, b, bnext, f, x)                                         \
    {                                                                          \
        __linkern_flipper_flip((x), (a), (b));                                 \
        (f)->stack[(f)->counter].first  = (a);                                 \
        (f)->stack[(f)->counter++].last = (b);                                 \
    }

#define UNFLIP(aprev, a, b, bnext, f, x)                                       \
    {                                                                          \
        __linkern_flipper_flip((x), (b), (a));                                 \
        (f)->counter--;                                                        \
    }

#define markedge_add(n1, n2, E) E->add_edges[n1 ^ n2] = 1
#define markedge_del(n1, n2, E) E->del_edges[n1 ^ n2] = 1
#define unmarkedge_add(n1, n2, E) E->add_edges[n1 ^ n2] = 0
#define unmarkedge_del(n1, n2, E) E->del_edges[n1 ^ n2] = 0
#define is_it_added(n1, n2, E) E->add_edges[n1 ^ n2]
#define is_it_deleted(n1, n2, E) E->del_edges[n1 ^ n2]

typedef struct edge
{
    int other;
    int weight;
} edge;

typedef struct edgelook edgelook;
typedef struct edgelook
{
    edgelook *next;
    int other;
    int diff;
    int over;
    int seq;
    int side;
    int mm;
} edgelook;

typedef struct flippair
{
    int firstprev;
    int first;
    int last;
    int lastnext;
} flippair;

typedef struct flipstack
{
    flippair *stack;
    int counter;
    int max;
} flipstack;

typedef struct lk_graph
{
    edge **goodlist;
    edge *edgespace;
    int *degree;
    int *weirdmark;
    int weirdmagic;
    int ncount;
} lk_graph;

typedef struct adddel
{
    char *add_edges;
    char *del_edges;
} adddel;

typedef struct active_node active_node;
typedef struct active_node
{
    int i;
    active_node *next;
} queue_item;

typedef struct active_queue
{
    char *active;
    active_node *active_queue;
    active_node *bottom_active_queue;
} active_queue;

static void
lin_kernighan(tsp_prob *tsp, tsp_heur_linkern_env *lk_env, lk_graph *G,
              adddel *E, active_queue *Q, __lk_flipper *F, double *val,
              tsp_sol *win_sol, flipstack *win, flipstack *fstack),
look_ahead_noback(lk_graph *G, solver_data *data, adddel *E, __lk_flipper *F,
                  int first, int last, int gain, edgelook *winner),
bigturn(lk_graph *G, int n, int tonext, active_queue *Q, __lk_flipper *F,
        solver_data *data),
first_kicker(lk_graph *G, solver_data *data, __lk_flipper *F, int *t1, int *t2),
find_close_four(lk_graph *G, solver_data *data, __lk_flipper *F, int *t1,
                int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8),
find_walk_four(lk_graph *G, solver_data *data, __lk_flipper *F, int *t1,
               int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8),
randcycle(int ncount, int *cyc), insertedge(lk_graph *G, int n1, int n2, int w),
initgraph(lk_graph *G), freegraph(lk_graph *G), free_flipstack(flipstack *f);

static void
adddel_init(adddel *E),
adddel_free(adddel *E), active_queue_init(active_queue *Q),
active_queue_free(active_queue *Q),
active_queue_add_node(active_queue *Q, int n);

static int
adddel_build(adddel *E, int ncount),
active_queue_build(active_queue *Q, int ncount),
active_queue_pop_head(active_queue *Q);

static int
buildgraph(lk_graph *G, int ncount, int ecount, int *elist, solver_data *data),
repeated_lin_kernighan(tsp_prob *tsp, tsp_heur_linkern_env *lk_env,
                       tsp_sol *sol, lk_graph *G),
weird_second_step(lk_graph *G, solver_data *data, adddel *E, active_queue *Q,
                  __lk_flipper *F, int gain, int t1, int t2, flipstack *fstack),
step(lk_graph *G, solver_data *data, adddel *E, active_queue *Q,
     __lk_flipper *F, int level, int gain, int *Gstar, int first, int last,
     flipstack *fstack),
step_noback(lk_graph *G, solver_data *data, adddel *E, active_queue *Q,
            __lk_flipper *F, int level, int gain, int *Gstar, int first,
            int last, flipstack *fstack),
find_geometric_four(lk_graph *G, solver_data *data, __lk_flipper *F, int *t1,
                    int *t2, int *t3, int *t4, int *t5, int *t6, int *t7,
                    int *t8),
random_four_swap(lk_graph *G, solver_data *data, active_queue *Q,
                 __lk_flipper *F, int *delta, int kicktype, flipstack *win,
                 flipstack *fstack),
init_flipstack(flipstack *f, int total, int single);

static double
improve_tour(lk_graph *G, solver_data *data, adddel *E, active_queue *Q,
             __lk_flipper *F, int start, flipstack *fstack);

static edgelook *
look_ahead(lk_graph *G, solver_data *data, adddel *E, __lk_flipper *F,
           int first, int last, int gain, int level),
*weird_look_ahead(lk_graph *G, solver_data *data, __lk_flipper *F, int gain,
                  int t1, int t2),
*weird_look_ahead2(lk_graph *G, solver_data *data, __lk_flipper *F, int gain,
                   int t2, int t3, int t4),
*weird_look_ahead3(lk_graph *G, solver_data *data, __lk_flipper *F, int gain,
                   int t2, int t3, int t6);

int
tsp_opt_heur_linkern(tsp_prob *tsp, tsp_heur_linkern_env *lk_env, tsp_sol *sol)
{
    int rval      = 0;
    tsp_sol *tsol = NULL;
    lk_graph G;
    tsp_heur_linkern_param *lk_param = lk_env->param;
    tsp_heur_linkern_stats *lk_stats = lk_env->stats;
    int kicktype_old                 = -1;
    int nkicks_old                   = -1;
    data_map *map                    = tsp->data->map;

    if (lk_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("linkern ...\n");
    }

    rval = stats_start(lk_stats->total);
    check_rval(rval, "stats_start failed", CLEANUP);

    initgraph(&G);

    if (map->img_n < 10 && lk_param->nkicks > 0)
    {
        if (lk_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("Less than 10 nodes, setting repeatcount to 0\n");
        }
        nkicks_old       = lk_param->nkicks;
        lk_param->nkicks = 0;
    }

    if (!data_is_norm_type(tsp->data, SOLVER_DATA_TYPE_EUCLIDEAN))
    {
        if (lk_param->kick_type == TSP_LINKERN_GEOMETRIC_KICK)
        {
            if (lk_env->verbosity >= SOLVER_VERBOSITY_INFO)
            {
                printf("Setting kick type to close\n");
            }
            kicktype_old        = lk_param->kick_type;
            lk_param->kick_type = TSP_LINKERN_CLOSE_KICK;
        }
    }

    rval = buildgraph(&G, map->img_n, map->kn_ecount, map->kn_elist, tsp->data);
    check_rval(rval, "buildgraph failed", CLEANUP);

    tsol = tsp_create_sol(tsp);
    check_null(tsol, "out of memory in linkern", CLEANUP);

    tsp_copy_sol(sol, tsol);

    if (lk_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("Starting Cycle: %.0f\n", tsol->val);
    }

    rval = repeated_lin_kernighan(tsp, lk_env, tsol, &G);
    check_rval(rval, "repeated_lin_kernighan failed", CLEANUP);

    rval = stats_stop(lk_stats->total, !rval);
    check_rval(rval, "stats_stop failed", CLEANUP);

    if (lk_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("Best cycle length: %.0f\n", tsol->val);
        printf("Lin-Kernighan Running Time: %.2ld\n",
               stats_get_current_time(lk_stats->total));
    }

    tsp_copy_sol(tsol, sol);

CLEANUP:

    if (kicktype_old != -1)
        lk_param->kick_type = kicktype_old;
    if (nkicks_old != -1)
        lk_param->nkicks = nkicks_old;

    tsp_free_sol(&tsol);
    freegraph(&G);

    return rval;
}

#define STALLCOUNT 100000000

static int
repeated_lin_kernighan(tsp_prob *tsp, tsp_heur_linkern_env *lk_env,
                       tsp_sol *sol, lk_graph *G)
{
    int rval  = 0;
    int round = 0;
    int quitcount, hit, delta;
    tsp_sol *win_sol = NULL;
    flipstack winstack, fstack;
    double t, best = sol->val;
    int ncount = G->ncount;
    adddel E;
    __lk_flipper F;
    active_queue Q;
    tsp_heur_linkern_param *lk_param = lk_env->param;
    tsp_heur_linkern_stats *lk_stats = lk_env->stats;
    int initial_count;

    initial_count = lk_stats->steps->count_active;

    rval = stats_start(lk_stats->steps);
    check_rval(rval, "stats_start failed", CLEANUP);

    active_queue_init(&Q);
    adddel_init(&E);

    rval = active_queue_build(&Q, ncount);
    check_rval(rval, "failed", CLEANUP);
    rval = adddel_build(&E, ncount);
    check_rval(rval, "failed", CLEANUP);

    hit  = 2 * (MAXDEPTH + 7 + KICK_MAXDEPTH);
    rval = init_flipstack(&fstack, hit, 0);
    check_rval(rval, "init_flipstack failed", CLEANUP);
    rval = init_flipstack(&winstack, 500 + ncount / 50, hit);
    check_rval(rval, "init_flipstack failed", CLEANUP);

    win_sol = tsp_create_sol(tsp);
    check_null(win_sol, "out of memory in repeated_lin_kernighan", CLEANUP);
    win_sol->cycle[0] = -1;

    quitcount = STALLCOUNT;
    if (quitcount > lk_param->nkicks)
        quitcount = lk_param->nkicks;

    __linkern_flipper_init(&F, sol);
    fstack.counter    = 0;
    winstack.counter  = 0;
    win_sol->cycle[0] = -1;

    int *tcyc = NULL;
    int i;

    tcyc = malloc(ncount * sizeof(int));
    check_null(tcyc, "out of memory in repeated_lin_kernighan", CLEANUP);
    /* init active_queue with random order */
    randcycle(ncount, tcyc);
    for (i = 0; i < ncount; i++)
    {
        active_queue_add_node(&Q, tcyc[i]);
    }
    free(tcyc);

    lin_kernighan(tsp, lk_env, G, &E, &Q, &F, &best, win_sol, &winstack,
                  &fstack);

    rval = stats_stop(lk_stats->steps, 1);
    check_rval(rval, "stats_stop failed", CLEANUP);

    winstack.counter  = 0;
    win_sol->cycle[0] = -1;

    if (lk_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        if (quitcount > 0)
        {
            printf("%4d Steps   Best: %.0f   %.2ld seconds\n",
                   lk_stats->steps->count_active - initial_count, best,
                   stats_get_total_time(lk_stats->steps) - 1);
        }
        else
        {
            printf("LK Cycle: %.0f\n", best);
        }
    }

    round = lk_stats->steps->count_active - initial_count;

    while (round < quitcount)
    {
        rval = stats_start(lk_stats->steps);
        check_rval(rval, "stats_start failed", CLEANUP);

        hit            = 0;
        fstack.counter = 0;

        rval = random_four_swap(G, tsp->data, &Q, &F, &delta,
                                lk_param->kick_type, &winstack, &fstack);
        check_rval(rval, "random_four_swap failed", CLEANUP);

        fstack.counter = 0;
        t              = best + delta;
        lin_kernighan(tsp, lk_env, G, &E, &Q, &F, &t, win_sol, &winstack,
                      &fstack);

        if (t <= best)
        {
            winstack.counter  = 0;
            win_sol->cycle[0] = -1;
            if (t < best)
            {
                best = t;
                quitcount =
                lk_stats->steps->count_active - initial_count + STALLCOUNT;
                if (quitcount > lk_param->nkicks)
                    quitcount = lk_param->nkicks;
                hit++;
            }
        }
        else
        {
            if (win_sol->cycle[0] == -1)
            {
                while (winstack.counter)
                {
                    winstack.counter--;
                    __linkern_flipper_flip(
                    &F, winstack.stack[winstack.counter].last,
                    winstack.stack[winstack.counter].first);
                }
            }
            else
            {
                __linkern_flipper_finish(&F);
                __linkern_flipper_init(&F, win_sol);
                while (winstack.counter)
                {
                    winstack.counter--;
                    __linkern_flipper_flip(
                    &F, winstack.stack[winstack.counter].last,
                    winstack.stack[winstack.counter].first);
                }
                win_sol->cycle[0] = -1;
            }
        }

        rval = stats_stop(lk_stats->steps, 1);
        check_rval(rval, "stats_stop failed", CLEANUP);

        round = lk_stats->steps->count_active - initial_count;
        if (lk_env->verbosity >= SOLVER_VERBOSITY_INFO &&
            (hit || (round % 1000 == 999)))
        {
            printf("%4d Steps   Best: %.0f   %.2ld seconds\n", round, best,
                   stats_get_current_time(lk_stats->steps));
        }

        if (lk_param->time_limit > 0.0 &&
            (stats_get_current_time(lk_stats->steps)) > lk_param->time_limit)
        {
            printf("STOP - timebound (%.2ld seconds)\n",
                   stats_get_current_time(lk_stats->steps));
            if (lk_env->verbosity >= SOLVER_VERBOSITY_INFO)
            {
                printf("STEPS: %d\n", round);
            }
            break;
        }
        if (lk_param->length_bound > 0.0 && best <= lk_param->length_bound)
        {
            printf("STOP - length bound reached (%.0f)\n",
                   lk_param->length_bound);
            if (lk_env->verbosity >= SOLVER_VERBOSITY_INFO)
            {
                printf("STEPS: %d\n", round);
            }
            break;
        }
    }

    if (lk_env->verbosity >= SOLVER_VERBOSITY_INFO && round > 0)
    {
        printf("%4d Total Steps.\n", round);
    }

    __linkern_flipper_cycle(&F, sol, tsp->data);
    __linkern_flipper_finish(&F);

    assert(best == sol->val);

CLEANUP:

    active_queue_free(&Q);
    adddel_free(&E);
    free_flipstack(&fstack);
    free_flipstack(&winstack);
    tsp_free_sol(&win_sol);
    return rval;
}

static void
lin_kernighan(tsp_prob *tsp, tsp_heur_linkern_env *lk_env, lk_graph *G,
              adddel *E, active_queue *Q, __lk_flipper *F, double *val,
              tsp_sol *win_sol, flipstack *win, flipstack *fstack)
{
    int start, i;
    double delta, totalwin = 0.0;

    while (1)
    {
        start = active_queue_pop_head(Q);
        if (start == -1)
            break;

        delta = improve_tour(G, tsp->data, E, Q, F, start, fstack);
        if (delta > 0.0)
        {
            totalwin += delta;
            if (win->counter < win->max)
            {
                for (i = 0; i < fstack->counter; i++)
                {
                    win->stack[win->counter].first = fstack->stack[i].first;
                    win->stack[win->counter].last  = fstack->stack[i].last;
                    win->stack[win->counter].firstprev =
                    fstack->stack[i].firstprev;
                    win->stack[win->counter].lastnext =
                    fstack->stack[i].lastnext;
                    win->counter++;
                }
            }
            else if (win_sol->cycle[0] == -1)
            {
                for (i = 0; i < fstack->counter; i++)
                {
                    win->stack[win->counter].first = fstack->stack[i].first;
                    win->stack[win->counter].last  = fstack->stack[i].last;
                    win->counter++;
                }
                __linkern_flipper_cycle(F, win_sol, tsp->data);
            }
            fstack->counter = 0;
        }
    }

    if (win_sol->cycle[0] == -1)
    {
        for (i = 0; i < fstack->counter; i++)
        {
            win->stack[win->counter].first = fstack->stack[i].first;
            win->stack[win->counter].last  = fstack->stack[i].last;
            win->counter++;
        }
    }
    (*val) -= totalwin;
}

static double
improve_tour(lk_graph *G, solver_data *data, adddel *E, active_queue *Q,
             __lk_flipper *F, int t1, flipstack *fstack)
{
    int t2 = __linkern_flipper_next(F, t1);
    int gain, Gstar = 0;

    gain = Edgelen(t1, t2, data);
    markedge_del(t1, t2, E);

    if (step(G, data, E, Q, F, 0, gain, &Gstar, t1, t2, fstack) == 0)
    {
        Gstar = weird_second_step(G, data, E, Q, F, gain, t1, t2, fstack);
    }
    unmarkedge_del(t1, t2, E);

    if (Gstar)
    {
        active_queue_add_node(Q, t1);
        active_queue_add_node(Q, t2);
    }
    return (double)Gstar;
}

static int
step(lk_graph *G, solver_data *data, adddel *E, active_queue *Q,
     __lk_flipper *F, int level, int gain, int *Gstar, int first, int last,
     flipstack *fstack)
{
    int val, this, newlast, hit = 0, oldG = gain;
    edgelook *list, *e, *efree, *next;

    if (level >= BACKTRACK)
    {
        return step_noback(G, data, E, Q, F, level, gain, Gstar, first, last,
                           fstack);
    }

    list = look_ahead(G, data, E, F, first, last, gain, level);
    for (e = list; e; e = e->next)
    {
        this    = e->other;
        newlast = e->over;

        gain = oldG - e->diff;
        val  = gain - Edgelen(newlast, first, data);
        if (val > *Gstar)
        {
            *Gstar = val;
            hit++;
        }

        FLIP(first, last, newlast, this, fstack, F);

        if (level < MAXDEPTH)
        {
            markedge_add(last, this, E);
            markedge_del(this, newlast, E);
            hit += step(G, data, E, Q, F, level + 1, gain, Gstar, first,
                        newlast, fstack);
            unmarkedge_add(last, this, E);
            unmarkedge_del(this, newlast, E);
        }

        if (!hit)
        {
            UNFLIP(first, last, newlast, this, fstack, F);
        }
        else
        {
            active_queue_add_node(Q, this);
            active_queue_add_node(Q, newlast);
            for (efree = list; efree; efree = next)
            {
                next = efree->next;
                free(efree);
            }
            return 1;
        }
    }

    for (efree = list; efree; efree = next)
    {
        next = efree->next;
        free(efree);
    }

    return 0;
}

static int
step_noback(lk_graph *G, solver_data *data, adddel *E, active_queue *Q,
            __lk_flipper *F, int level, int gain, int *Gstar, int first,
            int last, flipstack *fstack)
{
    edgelook e;

    look_ahead_noback(G, data, E, F, first, last, gain - *Gstar - level, &e);

    if (e.diff < SOLVER_MAXINT)
    {
        if (e.mm)
        {
            int hit      = 0;
            int this     = e.other;
            int newfirst = e.over;
            int val;

            gain -= e.diff;
            val = gain - Edgelen(newfirst, last, data);
            if (val > *Gstar)
            {
                *Gstar = val;
                hit++;
            }
            FLIP(this, newfirst, first, last, fstack, F);

            if (level < MAXDEPTH)
            {
                markedge_add(first, this, E);
                markedge_del(this, newfirst, E);
                hit += step_noback(G, data, E, Q, F, level + 1, gain, Gstar,
                                   newfirst, last, fstack);
                unmarkedge_add(first, this, E);
                unmarkedge_del(this, newfirst, E);
            }

            if (!hit)
            {
                UNFLIP(this, newfirst, first, last, fstack, F);
                return 0;
            }
            else
            {
                active_queue_add_node(Q, this);
                active_queue_add_node(Q, newfirst);
                return 1;
            }
        }
        else
        {
            int hit     = 0;
            int this    = e.other;
            int newlast = e.over;
            int val;

            gain -= e.diff;
            val = gain - Edgelen(newlast, first, data);
            if (val > *Gstar)
            {
                *Gstar = val;
                hit++;
            }

            FLIP(first, last, newlast, this, fstack, F);

            if (level < MAXDEPTH)
            {
                markedge_add(last, this, E);
                markedge_del(this, newlast, E);
                hit += step_noback(G, data, E, Q, F, level + 1, gain, Gstar,
                                   first, newlast, fstack);
                unmarkedge_add(last, this, E);
                unmarkedge_del(this, newlast, E);
            }

            if (!hit)
            {
                UNFLIP(first, last, newlast, this, fstack, F);
                return 0;
            }
            else
            {
                active_queue_add_node(Q, this);
                active_queue_add_node(Q, newlast);
                return 1;
            }
        }
    }
    else
    {
        return 0;
    }
}

#define G_MULT 1.5

static int
weird_second_step(lk_graph *G, solver_data *data, adddel *E, active_queue *Q,
                  __lk_flipper *F, int len_t1_t2, int t1, int t2,
                  flipstack *fstack)
{
    int t3, t4, t5, t6, t7, t8;
    int oldG, gain, tG, Gstar = 0, val, hit;
    int t3prev, t4next;
    edgelook *e, *f, *h, *list, *list2, *list3, *efree, *next;

    list = weird_look_ahead(G, data, F, len_t1_t2, t1, t2);
    for (h = list; h; h = h->next)
    {
        t3 = h->other;
        t4 = h->over;

        oldG = len_t1_t2 - h->diff;

        t3prev = __linkern_flipper_prev(F, t3);
        t4next = __linkern_flipper_next(F, t4);

        markedge_add(t2, t3, E);
        markedge_del(t3, t4, E);
        G->weirdmagic++;
        G->weirdmark[t1]     = G->weirdmagic;
        G->weirdmark[t2]     = G->weirdmagic;
        G->weirdmark[t3]     = G->weirdmagic;
        G->weirdmark[t4next] = G->weirdmagic;

        list2 = weird_look_ahead2(G, data, F, oldG, t2, t3, t4);
        for (e = list2; e; e = e->next)
        {
            t5 = e->other;
            t6 = e->over;

            markedge_add(t4, t5, E);
            if (e->seq)
            {
                if (!e->side)
                {
                    gain = oldG - e->diff;
                    val  = gain - Edgelen(t6, t1, data);
                    if (val > Gstar)
                        Gstar = val;
                    FLIP(t1, t2, t6, t5, fstack, F);
                    FLIP(t2, t5, t3, t4, fstack, F);

                    markedge_del(t5, t6, E);
                    hit =
                    step(G, data, E, Q, F, 2, gain, &Gstar, t1, t6, fstack);
                    unmarkedge_del(t5, t6, E);

                    if (!hit && Gstar)
                        hit = 1;

                    if (!hit)
                    {
                        UNFLIP(t2, t5, t3, t4, fstack, F);
                        UNFLIP(t1, t2, t6, t5, fstack, F);
                    }
                    else
                    {
                        unmarkedge_add(t2, t3, E);
                        unmarkedge_del(t3, t4, E);
                        unmarkedge_add(t4, t5, E);
                        active_queue_add_node(Q, t3);
                        active_queue_add_node(Q, t4);
                        active_queue_add_node(Q, t5);
                        active_queue_add_node(Q, t6);
                        for (e = list; e; e = next)
                        {
                            next = e->next;
                            free(e);
                        }
                        for (e = list2; e; e = next)
                        {
                            next = e->next;
                            free(e);
                        }
                        return Gstar;
                    }
                }
                else
                {
                    gain = oldG - e->diff;
                    val  = gain - Edgelen(t6, t1, data);
                    if (val > Gstar)
                        Gstar = val;
                    FLIP(t1, t2, t3, t4, fstack, F);
                    FLIP(t6, t5, t2, t4, fstack, F);
                    FLIP(t1, t3, t6, t2, fstack, F);

                    markedge_del(t5, t6, E);
                    hit =
                    step(G, data, E, Q, F, 2, gain, &Gstar, t1, t6, fstack);
                    unmarkedge_del(t5, t6, E);

                    if (!hit && Gstar)
                        hit = 1;

                    if (!hit)
                    {
                        UNFLIP(t1, t3, t6, t2, fstack, F);
                        UNFLIP(t6, t5, t2, t4, fstack, F);
                        UNFLIP(t1, t2, t3, t4, fstack, F);
                    }
                    else
                    {
                        unmarkedge_add(t2, t3, E);
                        unmarkedge_del(t3, t4, E);
                        unmarkedge_add(t4, t5, E);
                        active_queue_add_node(Q, t3);
                        active_queue_add_node(Q, t4);
                        active_queue_add_node(Q, t5);
                        active_queue_add_node(Q, t6);
                        for (e = list; e; e = next)
                        {
                            next = e->next;
                            free(e);
                        }
                        for (e = list2; e; e = next)
                        {
                            next = e->next;
                            free(e);
                        }
                        return Gstar;
                    }
                }
            }
            else
            {
                tG = oldG - e->diff;
                markedge_del(t5, t6, E);
                list3 = weird_look_ahead3(G, data, F, tG, t2, t3, t6);
                for (f = list3; f; f = f->next)
                {
                    t7   = f->other;
                    t8   = f->over;
                    gain = tG - f->diff;
                    if (!f->side)
                    {
                        val = gain - Edgelen(t8, t1, data);
                        if (val > Gstar)
                            Gstar = val;
                        FLIP(t1, t2, t8, t7, fstack, F);
                        FLIP(t2, t7, t3, t4, fstack, F);
                        FLIP(t7, t4, t6, t5, fstack, F);

                        markedge_add(t6, t7, E);
                        markedge_del(t7, t8, E);
                        hit =
                        step(G, data, E, Q, F, 3, gain, &Gstar, t1, t8, fstack);
                        unmarkedge_del(t6, t7, E);
                        unmarkedge_del(t7, t8, E);

                        if (!hit && Gstar)
                            hit = 1;

                        if (!hit)
                        {
                            UNFLIP(t7, t4, t6, t5, fstack, F);
                            UNFLIP(t2, t7, t3, t4, fstack, F);
                            UNFLIP(t1, t2, t8, t7, fstack, F);
                        }
                        else
                        {
                            unmarkedge_add(t2, t3, E);
                            unmarkedge_del(t3, t4, E);
                            unmarkedge_add(t4, t5, E);
                            unmarkedge_del(t5, t6, E);
                            active_queue_add_node(Q, t3);
                            active_queue_add_node(Q, t4);
                            active_queue_add_node(Q, t5);
                            active_queue_add_node(Q, t6);
                            active_queue_add_node(Q, t7);
                            active_queue_add_node(Q, t8);
                            for (e = list; e; e = next)
                            {
                                next = e->next;
                                free(e);
                            }
                            for (e = list2; e; e = next)
                            {
                                next = e->next;
                                free(e);
                            }
                            for (e = list3; e; e = next)
                            {
                                next = e->next;
                                free(e);
                            }
                            return Gstar;
                        }
                    }
                    else
                    {
                        val = gain - Edgelen(t8, t1, data);
                        if (val > Gstar)
                            Gstar = val;
                        FLIP(t1, t2, t6, t5, fstack, F);
                        FLIP(t1, t6, t8, t7, fstack, F);
                        FLIP(t3, t4, t2, t5, fstack, F);

                        markedge_add(t6, t7, E);
                        markedge_del(t7, t8, E);
                        hit =
                        step(G, data, E, Q, F, 3, gain, &Gstar, t1, t8, fstack);
                        unmarkedge_add(t6, t7, E);
                        unmarkedge_del(t7, t8, E);

                        if (!hit && Gstar)
                            hit = 1;

                        if (!hit)
                        {
                            UNFLIP(t3, t4, t2, t5, fstack, F);
                            UNFLIP(t1, t6, t8, t7, fstack, F);
                            UNFLIP(t1, t2, t6, t5, fstack, F);
                        }
                        else
                        {
                            unmarkedge_add(t2, t3, E);
                            unmarkedge_del(t3, t4, E);
                            unmarkedge_add(t4, t5, E);
                            unmarkedge_del(t5, t6, E);
                            active_queue_add_node(Q, t3);
                            active_queue_add_node(Q, t4);
                            active_queue_add_node(Q, t5);
                            active_queue_add_node(Q, t6);
                            active_queue_add_node(Q, t7);
                            active_queue_add_node(Q, t8);
                            for (e = list; e; e = next)
                            {
                                next = e->next;
                                free(e);
                            }
                            for (e = list2; e; e = next)
                            {
                                next = e->next;
                                free(e);
                            }
                            for (e = list3; e; e = next)
                            {
                                next = e->next;
                                free(e);
                            }
                            return Gstar;
                        }
                    }
                }
                for (efree = list3; efree; efree = next)
                {
                    next = efree->next;
                    free(efree);
                }
                unmarkedge_del(t5, t6, E);
            }
            unmarkedge_add(t4, t5, E);
        }
        for (efree = list2; efree; efree = next)
        {
            next = efree->next;
            free(efree);
        }
        unmarkedge_add(t2, t3, E);
        unmarkedge_del(t3, t4, E);
    }
    for (efree = list; efree; efree = next)
    {
        next = efree->next;
        free(efree);
    }
    return 0;
}

static edgelook *
look_ahead(lk_graph *G, solver_data *data, adddel *E, __lk_flipper *F,
           int first, int last, int gain, int level)
{
    edgelook *list = NULL, *el;
    int i, val;
    int this, prev;
    int lastnext = __linkern_flipper_next(F, last);
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, ahead = backtrack_count[level];
    edge **goodlist = G->goodlist;

    for (i = 0; i < ahead; i++)
    {
        value[i] = SOLVER_MAXINT;
    }
    value[ahead] = -SOLVER_MAXINT;

    for (i = 0; goodlist[last][i].weight <= gain; i++)
    {
        this = goodlist[last][i].other;
        if (!is_it_deleted(last, this, E) && this != first && this != lastnext)
        {
            prev = __linkern_flipper_prev(F, this);
            if (!is_it_added(this, prev, E))
            {
                val = goodlist[last][i].weight - Edgelen(this, prev, data);
                if (val < value[0])
                {
                    for (k = 0; value[k + 1] > val; k++)
                    {
                        value[k] = value[k + 1];
                        other[k] = other[k + 1];
                        save[k]  = save[k + 1];
                    }
                    value[k] = val;
                    other[k] = this;
                    save[k]  = prev;
                }
            }
        }
    }

    for (i = 0; i < ahead; i++)
    {
        if (value[i] < SOLVER_MAXINT)
        {
            el        = malloc(sizeof(edgelook));
            el->diff  = value[i];
            el->other = other[i];
            el->over  = save[i];
            el->next  = list;
            list      = el;
        }
    }

    return list;
}

static void
look_ahead_noback(lk_graph *G, solver_data *data, adddel *E, __lk_flipper *F,
                  int first, int last, int gain, edgelook *winner)
{
    int val;
    int this, prev;
    int lastnext = __linkern_flipper_next(F, last);
    int i;
    edge **goodlist = G->goodlist;

    winner->diff = SOLVER_MAXINT;
    for (i = 0; goodlist[last][i].weight < gain; i++)
    {
        this = goodlist[last][i].other;

        if (!is_it_deleted(last, this, E) && this != first && this != lastnext)
        {
            prev = __linkern_flipper_prev(F, this);
            if (!is_it_added(this, prev, E))
            {
                val = goodlist[last][i].weight - Edgelen(this, prev, data);
                if (val < winner->diff)
                {
                    winner->diff  = val;
                    winner->other = this;
                    winner->over  = prev;
                    winner->mm    = 0;
                }
            }
        }
    }
    int firstprev = __linkern_flipper_prev(F, first);
    int next;

    for (i = 0; goodlist[first][i].weight < gain; i++)
    {
        this = goodlist[first][i].other;
        if (!is_it_deleted(first, this, E) && this != last && this != firstprev)
        {
            next = __linkern_flipper_next(F, this);
            if (!is_it_added(this, next, E))
            {
                val = goodlist[first][i].weight - Edgelen(this, next, data);
                if (val < winner->diff)
                {
                    winner->diff  = val;
                    winner->other = this;
                    winner->over  = next;
                    winner->mm    = 1;
                }
            }
        }
    }
}

static edgelook *
weird_look_ahead(lk_graph *G, solver_data *data, __lk_flipper *F, int gain,
                 int t1, int t2)
{
    edgelook *list, *el;
    int i, this, next;
    int other[MAX_BACK], save[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val, ahead;
    edge **goodlist = G->goodlist;

    list  = NULL;
    ahead = weird_backtrack_count[0];
    for (i = 0; i < ahead; i++) value[i] = SOLVER_MAXINT;
    value[ahead] = -SOLVER_MAXINT;

    for (i = 0; goodlist[t2][i].weight <= gain; i++)
    {
        this = goodlist[t2][i].other;
        if (this != t1)
        {
            next = __linkern_flipper_next(F, this);
            val  = goodlist[t2][i].weight - Edgelen(this, next, data);
            if (val < value[0])
            {
                for (k = 0; value[k + 1] > val; k++)
                {
                    value[k] = value[k + 1];
                    other[k] = other[k + 1];
                    save[k]  = save[k + 1];
                }
                value[k] = val;
                other[k] = this;
                save[k]  = next;
            }
        }
    }
    for (i = 0; i < ahead; i++)
    {
        if (value[i] < SOLVER_MAXINT)
        {
            el        = malloc(sizeof(edgelook));
            el->diff  = value[i];
            el->other = other[i];
            el->over  = save[i];
            el->next  = list;
            list      = el;
        }
    }
    return list;
}

static edgelook *
weird_look_ahead2(lk_graph *G, solver_data *data, __lk_flipper *F, int gain,
                  int t2, int t3, int t4)
{
    edgelook *list = NULL;
    edgelook *el;
    int i, t5, t6;
    int other[MAX_BACK], save[MAX_BACK], seq[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val;
    int ahead       = weird_backtrack_count[1];
    edge **goodlist = G->goodlist;
    int *weirdmark  = G->weirdmark;
    int weirdmagic  = G->weirdmagic;

    for (i = 0; i < ahead; i++) value[i] = SOLVER_MAXINT;
    value[ahead] = -SOLVER_MAXINT;

    for (i = 0; goodlist[t4][i].weight <= gain; i++)
    {
        t5 = goodlist[t4][i].other;
        if (weirdmark[t5] != weirdmagic)
        {
            if (__linkern_flipper_sequence(F, t2, t5, t3))
            {
                t6  = __linkern_flipper_prev(F, t5);
                val = goodlist[t4][i].weight - Edgelen(t5, t6, data);
                if (val < value[0])
                {
                    for (k = 0; value[k + 1] > val; k++)
                    {
                        value[k] = value[k + 1];
                        other[k] = other[k + 1];
                        save[k]  = save[k + 1];
                        seq[k]   = seq[k + 1];
                        side[k]  = side[k + 1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k]  = t6;
                    seq[k]   = 1;
                    side[k]  = 0;
                }
                t6  = __linkern_flipper_next(F, t5);
                val = goodlist[t4][i].weight - Edgelen(t5, t6, data);
                if (val < value[0])
                {
                    for (k = 0; value[k + 1] > val; k++)
                    {
                        value[k] = value[k + 1];
                        other[k] = other[k + 1];
                        save[k]  = save[k + 1];
                        seq[k]   = seq[k + 1];
                        side[k]  = side[k + 1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k]  = t6;
                    seq[k]   = 1;
                    side[k]  = 1;
                }
            }
            else
            {
                t6  = __linkern_flipper_prev(F, t5);
                val = goodlist[t4][i].weight - Edgelen(t5, t6, data);
                if (val < value[0])
                {
                    for (k = 0; value[k + 1] > val; k++)
                    {
                        value[k] = value[k + 1];
                        other[k] = other[k + 1];
                        save[k]  = save[k + 1];
                        seq[k]   = seq[k + 1];
                        side[k]  = side[k + 1];
                    }
                    value[k] = val;
                    other[k] = t5;
                    save[k]  = t6;
                    seq[k]   = 0;
                    side[k]  = 0;
                }
            }
        }
    }

    for (i = 0; i < ahead; i++)
    {
        if (value[i] < SOLVER_MAXINT)
        {
            el        = malloc(sizeof(edgelook));
            el->diff  = value[i];
            el->other = other[i];
            el->over  = save[i];
            el->seq   = seq[i];
            el->side  = side[i];
            el->next  = list;
            list      = el;
        }
    }
    return list;
}

static edgelook *
weird_look_ahead3(lk_graph *G, solver_data *data, __lk_flipper *F, int gain,
                  int t2, int t3, int t6)
{
    edgelook *list = NULL;
    edgelook *el;
    int i, t7, t8;
    int other[MAX_BACK], save[MAX_BACK], side[MAX_BACK];
    int value[MAX_BACK + 1];
    int k, val;
    int ahead       = weird_backtrack_count[2];
    edge **goodlist = G->goodlist;
    int *weirdmark  = G->weirdmark;
    int weirdmagic  = G->weirdmagic;

    for (i = 0; i < ahead; i++) value[i] = SOLVER_MAXINT;
    value[ahead] = -SOLVER_MAXINT;

    for (i = 0; goodlist[t6][i].weight <= gain; i++)
    {
        t7 = goodlist[t6][i].other; /* Need t7 != t2, t3, t2next, t3prev */
        if (weirdmark[t7] != weirdmagic &&
            __linkern_flipper_sequence(F, t2, t7, t3))
        {
            t8  = __linkern_flipper_prev(F, t7);
            val = goodlist[t6][i].weight - Edgelen(t7, t8, data);
            if (val < value[0])
            {
                for (k = 0; value[k + 1] > val; k++)
                {
                    value[k] = value[k + 1];
                    other[k] = other[k + 1];
                    save[k]  = save[k + 1];
                    side[k]  = side[k + 1];
                }
                value[k] = val;
                other[k] = t7;
                save[k]  = t8;
                side[k]  = 0;
            }
            t8  = __linkern_flipper_next(F, t7);
            val = goodlist[t6][i].weight - Edgelen(t7, t8, data);
            if (val < value[0])
            {
                for (k = 0; value[k + 1] > val; k++)
                {
                    value[k] = value[k + 1];
                    other[k] = other[k + 1];
                    save[k]  = save[k + 1];
                    side[k]  = side[k + 1];
                }
                value[k] = val;
                other[k] = t7;
                save[k]  = t8;
                side[k]  = 1;
            }
        }
    }

    for (i = 0; i < ahead; i++)
    {
        if (value[i] < SOLVER_MAXINT)
        {
            el        = malloc(sizeof(edgelook));
            el->diff  = value[i];
            el->other = other[i];
            el->over  = save[i];
            el->side  = side[i];
            el->next  = list;
            list      = el;
        }
    }
    return list;
}

static int
random_four_swap(lk_graph *G, solver_data *data, active_queue *Q,
                 __lk_flipper *F, int *delta, int kicktype, flipstack *win,
                 flipstack *fstack)
{
    int rval = 0;
    int t1, t2, t3, t4, t5, t6, t7, t8, temp;

    switch (kicktype)
    {
    case TSP_LINKERN_WALK_KICK:
        find_walk_four(G, data, F, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8);
        break;
    case TSP_LINKERN_CLOSE_KICK:
        find_close_four(G, data, F, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8);
        break;
    case TSP_LINKERN_GEOMETRIC_KICK:
        rval =
        find_geometric_four(G, data, F, &t1, &t2, &t3, &t4, &t5, &t6, &t7, &t8);
        if (rval)
        {
            fprintf(stderr, "find_geometric_four failed\n");
            return 1;
        }
        break;
    default:
        fprintf(stderr, "unknown kick type %d\n", kicktype);
        return 1;
    }

    if (!__linkern_flipper_sequence(F, t1, t3, t5))
    {
        SWAP(t3, t5, temp);
        SWAP(t4, t6, temp);
    }
    if (!__linkern_flipper_sequence(F, t1, t5, t7))
    {
        SWAP(t5, t7, temp);
        SWAP(t6, t8, temp);
        if (!__linkern_flipper_sequence(F, t1, t3, t5))
        {
            SWAP(t3, t5, temp);
            SWAP(t4, t6, temp);
        }
    }
    FLIP(t1, t2, t5, t6, fstack, F);
    FLIP(t4, t3, t7, t8, fstack, F);
    FLIP(t1, t5, t6, t2, fstack, F);

    if (win->counter < win->max)
    {
        win->stack[win->counter].first = t2;
        win->stack[win->counter].last  = t5;
        win->counter++;
    }
    if (win->counter < win->max)
    {
        win->stack[win->counter].first = t3;
        win->stack[win->counter].last  = t7;
        win->counter++;
    }
    if (win->counter < win->max)
    {
        win->stack[win->counter].first = t5;
        win->stack[win->counter].last  = t6;
        win->counter++;
    }

    bigturn(G, t1, 0, Q, F, data);
    bigturn(G, t2, 1, Q, F, data);
    bigturn(G, t3, 0, Q, F, data);
    bigturn(G, t4, 1, Q, F, data);
    bigturn(G, t5, 0, Q, F, data);
    bigturn(G, t6, 1, Q, F, data);
    bigturn(G, t7, 0, Q, F, data);
    bigturn(G, t8, 1, Q, F, data);

    *delta = Edgelen(t1, t6, data) + Edgelen(t2, t5, data) +
             Edgelen(t3, t8, data) + Edgelen(t4, t7, data) -
             Edgelen(t1, t2, data) - Edgelen(t3, t4, data) -
             Edgelen(t5, t6, data) - Edgelen(t7, t8, data);

CLEANUP:
    return 0;
}

#define HUNT_PORTION_LONG 0.001

static void
first_kicker(lk_graph *G, solver_data *data, __lk_flipper *F, int *t1, int *t2)
{
    int longcount = (int)((double)G->ncount * HUNT_PORTION_LONG) + 10;
    int i, best, try1, len, next, prev, nextl, prevl;
    int ncount      = G->ncount;
    edge **goodlist = G->goodlist;

    try1  = rand() % ncount;
    next  = __linkern_flipper_next(F, try1);
    prev  = __linkern_flipper_prev(F, try1);
    nextl = Edgelen(try1, next, data);
    prevl = Edgelen(try1, prev, data);
    if (nextl >= prevl)
    {
        *t1  = try1;
        *t2  = next;
        best = nextl - goodlist[*t1][0].weight;
    }
    else
    {
        *t1  = prev;
        *t2  = try1;
        best = prevl - goodlist[*t1][0].weight;
    }

    for (i = 0; i < longcount; i++)
    {
        try1  = rand() % ncount;
        next  = __linkern_flipper_next(F, try1);
        prev  = __linkern_flipper_prev(F, try1);
        nextl = Edgelen(try1, next, data);
        prevl = Edgelen(try1, prev, data);
        if (nextl >= prevl)
        {
            len = nextl - goodlist[try1][0].weight;
            if (len > best)
            {
                *t1 = try1;
                *t2 = next;
            }
        }
        else
        {
            len = prevl - goodlist[try1][0].weight;
            if (len > best)
            {
                *t1 = prev;
                *t2 = try1;
            }
        }
    }
}

#define HUNT_PORTION 0.03
#define RAND_TRYS 6 /* To find the 3 other edges */

static void
find_close_four(lk_graph *G, solver_data *data, __lk_flipper *F, int *t1,
                int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8)
{
    int s1, s2, s3, s4, s5, s6, s7, s8;
    int i, k, try1, trydist;
    int count = (int)((double)G->ncount * HUNT_PORTION) + 1 + RAND_TRYS;
    int trials[RAND_TRYS + 1];
    int tdist[RAND_TRYS + 1];

    first_kicker(G, data, F, &s1, &s2);

TRYAGAIN:

    for (k = 0; k < RAND_TRYS; k++) tdist[k] = SOLVER_MAXINT;
    tdist[RAND_TRYS] = -SOLVER_MAXINT;
    for (i = 0; i < count; i++)
    {
        try1    = rand() % G->ncount;
        trydist = Edgelen(try1, s1, data);
        if (trydist < tdist[0])
        {
            for (k = 0; tdist[k + 1] > trydist; k++)
            {
                tdist[k]  = tdist[k + 1];
                trials[k] = trials[k + 1];
            }
            tdist[k]  = trydist;
            trials[k] = try1;
        }
    }

    k = RAND_TRYS - 1;
    do
    {
        if (k < 0)
            goto TRYAGAIN;
        s3 = trials[k--];
        s4 = __linkern_flipper_next(F, s3);
    } while (s3 == s1 || s3 == s2 || s4 == s1);

    do
    {
        if (k < 0)
            goto TRYAGAIN;
        s5 = trials[k--];
        s6 = __linkern_flipper_next(F, s5);
    } while (s5 == s1 || s5 == s2 || s5 == s3 || s5 == s4 || s6 == s1 ||
             s6 == s3);

    do
    {
        if (k < 0)
            goto TRYAGAIN;
        s7 = trials[k--];
        s8 = __linkern_flipper_next(F, s7);
    } while (s7 == s1 || s7 == s2 || s7 == s3 || s7 == s4 || s7 == s5 ||
             s7 == s6 || s8 == s1 || s8 == s3 || s8 == s5);

    *t1 = s1;
    *t2 = s2;
    *t3 = s3;
    *t4 = s4;
    *t5 = s5;
    *t6 = s6;
    *t7 = s7;
    *t8 = s8;
}

#define GEO_FACTOR 50
#define GEO_MAX 250

static int
find_geometric_four(lk_graph *G, solver_data *data, __lk_flipper *F, int *t1,
                    int *t2, int *t3, int *t4, int *t5, int *t6, int *t7,
                    int *t8)
{
    int rval   = 0;
    int *neigh = NULL;
    int temp, i, k, s1, s2, s3, s4, s5, s6, s7, s8;
    int trys;

    first_kicker(G, data, F, &s1, &s2);
    trys = (G->ncount / GEO_FACTOR) + 25;
    if (trys > GEO_MAX)
        trys = GEO_MAX;
    if (trys > G->ncount - 1)
        trys = G->ncount - 1;

    rval = data_get_node_k_nearest_euclidean(data, s1, trys, &neigh);
    check_rval(rval, " failed", CLEANUP);

    for (i = trys; i > trys - 9; i--)
    {
        k = rand() % i;
        SWAP(neigh[i - 1], neigh[k], temp);
    }

    k = trys - 1;
    do
    {
        s3 = neigh[k--];
        s4 = __linkern_flipper_next(F, s3);
    } while (s3 == s2 || s4 == s1);

    do
    {
        s5 = neigh[k--];
        s6 = __linkern_flipper_next(F, s5);
    } while (s5 == s2 || s5 == s4 || s6 == s1 || s6 == s3);

    do
    {
        s7 = neigh[k--];
        s8 = __linkern_flipper_next(F, s7);
    } while (s7 == s2 || s7 == s4 || s7 == s6 || s8 == s1 || s8 == s3 ||
             s8 == s5);

    *t1 = s1;
    *t2 = s2;
    *t3 = s3;
    *t4 = s4;
    *t5 = s5;
    *t6 = s6;
    *t7 = s7;
    *t8 = s8;

CLEANUP:
    if (neigh)
        free(neigh);

    return rval;
}

#define WALK_STEPS 50

static void
find_walk_four(lk_graph *G, solver_data *data, __lk_flipper *F, int *t1,
               int *t2, int *t3, int *t4, int *t5, int *t6, int *t7, int *t8)
{
    int s1, s2, s3, s4, s5, s6, s7, s8;
    int old, n, i, j;

    first_kicker(G, data, F, &s1, &s2);

    do
    {
        old = -1;
        n   = s2;

        for (i = 0; i < WALK_STEPS; i++)
        {
            j = rand() % (G->degree[n]);
            if (old != G->goodlist[n][j].other)
            {
                old = n;
                n   = G->goodlist[n][j].other;
            }
        }
        s3 = n;
        s4 = __linkern_flipper_next(F, s3);

        n = s4;
        for (i = 0; i < WALK_STEPS; i++)
        {
            j = rand() % (G->degree[n]);
            if (old != G->goodlist[n][j].other)
            {
                old = n;
                n   = G->goodlist[n][j].other;
            }
        }
        s5 = n;
        s6 = __linkern_flipper_next(F, s5);

        n = s6;
        for (i = 0; i < WALK_STEPS; i++)
        {
            j = rand() % (G->degree[n]);
            if (old != G->goodlist[n][j].other)
            {
                old = n;
                n   = G->goodlist[n][j].other;
            }
        }
        s7 = n;
        s8 = __linkern_flipper_next(F, s7);
    } while (s1 == s3 || s1 == s4 || s1 == s5 || s1 == s6 || s1 == s7 ||
             s1 == s8 || s2 == s3 || s2 == s4 || s2 == s5 || s2 == s6 ||
             s2 == s7 || s2 == s8 || s3 == s5 || s3 == s6 || s3 == s7 ||
             s3 == s8 || s4 == s5 || s4 == s6 || s4 == s7 || s4 == s8 ||
             s5 == s7 || s5 == s8 || s6 == s7 || s6 == s8);

    *t1 = s1;
    *t2 = s2;
    *t3 = s3;
    *t4 = s4;
    *t5 = s5;
    *t6 = s6;
    *t7 = s7;
    *t8 = s8;
}

static void
bigturn(lk_graph *G, int n, int tonext, active_queue *Q, __lk_flipper *F,
        solver_data *data)
{
    int i, k;

    active_queue_add_node(Q, n);
    if (tonext)
    {
        for (i = 0, k = n; i < MARK_LEVEL; i++)
        {
            k = __linkern_flipper_next(F, k);
            active_queue_add_node(Q, k);
        }
    }
    else
    {
        for (i = 0, k = n; i < MARK_LEVEL; i++)
        {
            k = __linkern_flipper_prev(F, k);
            active_queue_add_node(Q, k);
        }
    }

    for (i = 0; i < G->degree[n]; i++)
    {
        active_queue_add_node(Q, G->goodlist[n][i].other);
    }
}

static void
randcycle(int ncount, int *cyc)
{
    int i, k, temp;

    for (i = 0; i < ncount; i++) cyc[i] = i;
    for (i = ncount; i > 1; i--)
    {
        k = rand() % i;
        SWAP(cyc[i - 1], cyc[k], temp);
    }
}

static void
initgraph(lk_graph *G)
{
    G->goodlist   = NULL;
    G->edgespace  = NULL;
    G->degree     = NULL;
    G->weirdmark  = NULL;
    G->weirdmagic = 0;
    G->ncount     = 0;
}

static void
freegraph(lk_graph *G)
{
    if (G)
    {
        free(G->goodlist);
        free(G->edgespace);
        free(G->degree);
        free(G->weirdmark);
        G->weirdmagic = 0;
        G->ncount     = 0;
    }
}

static int
buildgraph(lk_graph *G, int ncount, int ecount, int *elist, solver_data *data)
{
    int rval = 0;
    int n1, n2, w, i;
    edge *p;

    G->goodlist = malloc(ncount * sizeof(edge *));
    check_null(G->goodlist, "out of memory in buildgraph", CLEANUP);
    G->degree = malloc(ncount * sizeof(int));
    check_null(G->degree, "out of memory in buildgraph", CLEANUP);
    G->weirdmark = malloc(ncount * sizeof(int));
    check_null(G->weirdmark, "out of memory in buildgraph", CLEANUP);
    G->edgespace = malloc(((2 * ecount) + ncount) * sizeof(edge));
    check_null(G->edgespace, "out of memory in buildgraph", CLEANUP);

    for (i = 0; i < ncount; i++)
    {
        G->degree[i]    = 1;
        G->weirdmark[i] = 0;
    }
    for (i = ecount - 1; i >= 0; i--)
    {
        G->degree[elist[2 * i]]++;
        G->degree[elist[(2 * i) + 1]]++;
    }

    for (i = 0, p = G->edgespace; i < ncount; i++)
    {
        G->goodlist[i] = p;
        p += (G->degree[i]);
        G->goodlist[i][G->degree[i] - 1].weight = SOLVER_MAXINT;
        G->degree[i]                            = 0;
    }

    assert(ecount);
    for (i = ecount - 1; i >= 0; i--)
    {
        n1 = elist[2 * i];
        n2 = elist[(2 * i) + 1];

        w = Edgelen(n1, n2, data);
        insertedge(G, n1, n2, w);
        insertedge(G, n2, n1, w);
    }
    G->ncount     = ncount;
    G->weirdmagic = 0;

CLEANUP:

    if (rval)
        freegraph(G);
    return rval;
}

static void
insertedge(lk_graph *G, int n1, int n2, int w)
{
    int i;
    edge *e = G->goodlist[n1];

    for (i = G->degree[n1] - 1; i >= 0 && e[i].weight >= w; i--)
    {
        e[i + 1].weight = e[i].weight;
        e[i + 1].other  = e[i].other;
    }
    e[i + 1].weight = w;
    e[i + 1].other  = n2;
    G->degree[n1]++;
}

static int
init_flipstack(flipstack *f, int total, int single)
{
    int rval   = 0;
    f->counter = 0;
    f->max     = 0;
    f->stack   = NULL;

    f->stack = malloc((total + single) * sizeof(flippair));
    check_null(f->stack, "out of memory in init_flipstack", CLEANUP);
    f->max = total;

CLEANUP:
    return rval;
}

static void
free_flipstack(flipstack *f)
{
    f->counter = 0;
    f->max     = 0;
    free(f->stack);
}

static void
adddel_init(adddel *E)
{
    E->add_edges = NULL;
    E->del_edges = NULL;
}

static void
adddel_free(adddel *E)
{
    if (E)
    {
        free(E->add_edges);
        free(E->del_edges);
    }
}

static int
adddel_build(adddel *E, int ncount)
{
    int rval = 0;
    int i, M;

    i = 0;
    while ((1 << i) < ncount) i++;
    M = (1 << i);

    E->add_edges = malloc(M * sizeof(char));
    check_null(E->add_edges, "out of memory ", CLEANUP);
    E->del_edges = malloc(M * sizeof(char));
    check_null(E->del_edges, "out of memory ", CLEANUP);

    for (i = 0; i < M; i++)
    {
        E->add_edges[i] = 0;
        E->del_edges[i] = 0;
    }

CLEANUP:

    if (rval)
    {
        adddel_free(E);
    }
    return rval;
}

static void
active_queue_init(active_queue *Q)
{
    Q->active              = NULL;
    Q->active_queue        = NULL;
    Q->bottom_active_queue = NULL;
}

static void
active_queue_free(active_queue *Q)
{
    active_node *node, *next;
    if (Q)
    {
        free(Q->active);
        for (node = Q->active_queue; node; node = next)
        {
            next = node->next;
            free(node);
        }
        Q->active_queue        = NULL;
        Q->bottom_active_queue = NULL;
    }
}

static int
active_queue_build(active_queue *Q, int ncount)
{
    int rval = 0;
    int i;

    active_queue_init(Q);

    Q->active = malloc(ncount * sizeof(char));
    check_null(Q->active, "out of memory", CLEANUP);

    for (i = 0; i < ncount; i++) Q->active[i] = 0;

CLEANUP:

    if (rval)
    {
        active_queue_free(Q);
    }
    return rval;
}

static void
active_queue_add_node(active_queue *Q, int n)
{
    active_node *node;

    if (Q->active[n] == 0)
    {
        Q->active[n] = 1;
        node         = malloc(sizeof(active_node));
        node->i      = n;
        node->next   = NULL;
        if (Q->bottom_active_queue)
        {
            Q->bottom_active_queue->next = node;
        }
        else
        {
            Q->active_queue = node;
        }
        Q->bottom_active_queue = node;
    }
}

static int
active_queue_pop_head(active_queue *Q)
{
    active_node *node;
    int n = -1;

    if (Q->active_queue != NULL)
    {
        node            = Q->active_queue;
        n               = node->i;
        Q->active_queue = node->next;
        if (node == Q->bottom_active_queue)
        {
            Q->bottom_active_queue = NULL;
        }
        free(node);
        Q->active[n] = 0;
    }
    return n;
}
