#include "ip/exact/bac/bac.h"
#include "cp/init/init.h"
#include "cp/heur/heur.h"
#include "cp/exact/bac/bac.h"

static void
randcycle(int ncount, int *cyc);

int
cp_get_xheur_ea__(void *cp_, void *bac_env_, void *sol_, int *nadded)
{
    int rval;
    cp_prob *cp               = (cp_prob *)cp_;
    cp_exact_bac_env *bac_env = (cp_exact_bac_env *)bac_env_;
    cp_sol *sol               = (cp_sol *)sol_;
    lp_prob *lp               = bac_env->ip->lp;
    cp_sol *tmpsol            = NULL;
    cp_cut *cuts              = NULL;

    tmpsol = cp_create_sol(cp);
    check_null(tmpsol, " out of memory", CLEANUP);

    rval = cp_get_xheur_ea(cp, bac_env, tmpsol);
    check_rval(rval, "failed", CLEANUP);

    if (bac_env->ip->lowerboundN < tmpsol->val - SOLVER_ZEROPLUS)
        bac_env->ip->lowerboundN = tmpsol->val;
    if (cp->ip->lowerboundG < tmpsol->val - SOLVER_ZEROPLUS)
    {
        cp->ip->lowerboundG = tmpsol->val;
        cp_copy_sol(tmpsol, sol);
    }

    cuts = cp_conv_cut_sol2connect(tmpsol);
    rval = cp_add_lp_cuts(cp, bac_env, &cuts, nadded);

    if (lp->status == SOLVER_LP_SUCCESS &&
        (bac_env->ip->upperboundN <
         bac_env->ip->lowerboundN + bac_env->ip->param->pruning_tol))
    {
        int edge_added;
        rval = cp_add_badvars(cp, bac_env, &edge_added);
    }

CLEANUP:
    cp_free_sol(&tmpsol);
    return rval;
}

int
cp_get_xheur_ea(cp_prob *cp, cp_exact_bac_env *bac_env, cp_sol *sol)
{
    int rval = 0;
    int i;

    op_pop *pop;
    lp_prob *lp           = bac_env->ip->lp;
    cp_heur_env *heur_env = bac_env->heur;

    graph_vertex *v;

    rval = stats_start(bac_env->ip->stats->xheur);

    int j, k, ns;
    if (bac_env->param->xheur_vph_meta)
        pop = cp_create_pop(cp, 10);
    else
        pop = cp_create_pop(cp, 1);

    for (i = 0; i < pop->size; i++)
    {
        cp_sol *tmpsol = pop->sol[i];
        cp_erase_sol(tmpsol);

        ns = 0;
        if (lp->sol->graph->na > 3)
        {

            do
            {
                for (k = 0; k < cp->n; k++)
                {
                    v = lp->sol->graph->v[k];
                    if (k == cp->data->from && !tmpsol->selected[k])
                    {
                        tmpsol->selected[cp->data->from] = 1;
                        tmpsol->cycle[ns++]              = cp->data->from;
                    }
                    else if (k == cp->data->to && !tmpsol->selected[k])
                    {
                        tmpsol->selected[cp->data->to] = 1;
                        tmpsol->cycle[ns++]            = cp->data->to;
                    }
                    else
                    {
                        if (SOLVER_ONEMINUS < v->y && !tmpsol->selected[k])
                        {
                            tmpsol->selected[k] = 1;
                            tmpsol->cycle[ns++] = k;
                        }
                        else if (SOLVER_ZEROPLUS < v->y && !tmpsol->selected[k])
                        {
                            if (rng_bernoulli(v->y))
                            {
                                tmpsol->selected[k] = 1;
                                tmpsol->cycle[ns++] = k;
                            }
                            else
                            {
                                tmpsol->selected[k] = 0;
                            }
                        }
                        else
                        {
                            tmpsol->selected[k] = 0;
                        }
                    }
                }
            } while (ns <= 3);
        }
        else
        {
            ;
        }

        tmpsol->ns  = ns;
        tmpsol->val = 0.0;

        for (j = 1; j < ns; j++)
        {
            tmpsol->cod_fr[tmpsol->cycle[j - 1]] = tmpsol->cycle[j];
            tmpsol->cod_bk[tmpsol->cycle[j]]     = tmpsol->cycle[j - 1];
            tmpsol->cap +=
            data_get_norm(cp->data, tmpsol->cycle[j - 1], tmpsol->cycle[j]);
            tmpsol->val += cp->data->obj_node[tmpsol->cycle[j - 1]];
        }
        tmpsol->cod_bk[tmpsol->cycle[0]]      = tmpsol->cycle[ns - 1];
        tmpsol->cod_fr[tmpsol->cycle[ns - 1]] = tmpsol->cycle[0];
        tmpsol->cap +=
        data_get_norm(cp->data, tmpsol->cycle[0], tmpsol->cycle[ns - 1]);
        tmpsol->val += cp->data->obj_node[tmpsol->cycle[ns - 1]];

        cp_init_sol(cp, bac_env->heur, tmpsol);
        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf(" %d: nvis: %d, length %.0f, value %.0f\n", i, tmpsol->ns,
                   tmpsol->cap, tmpsol->val);
        }
    }

    cp_update_pop(pop);
    cp_copy_sol(pop->sol[0], sol);

    if (bac_env->param->xheur_vph_meta)
    {
        cp_opt_heur_ea(cp, heur_env, pop, sol);
    }

    op_free_pop(&pop);

    rval = stats_stop(bac_env->ip->stats->xheur, 1);

    return rval;
}

static int
sort_edges(const void *xx, const void *yy)
{
    graph_arc *x = *(graph_arc **)xx, *y = *(graph_arc **)yy;

    if (x->x < y->x)
        return 1;
    if (x->x > y->x)
        return -1;

    return 0;
}

static void
randcycle(int ncount, int *cyc)
{
    int i, k, temp;
    for (i = 0; i < ncount; i++) cyc[i] = i;

    for (i = ncount; i > 2; i--)
    {
        k = rand() % i;
        SWAP(cyc[i - 1], cyc[k], temp);
    }
    return;
}

static int
call_random_tour(cp_prob *cp, solver_graph *graph, int nextremes,
                 graph_vertex **extremes)
{
    int npaths = ceil(nextremes / 2.0);
    int *shortcyc;
    int i, current;
    int zero_cnt = 0;
    graph_arc *arc;

    shortcyc = malloc(npaths * sizeof(int));
    for (i = 0; i < npaths; i++) shortcyc[i] = i;
    randcycle(npaths, shortcyc);

    if (!(graph->v[0]->deg))
        zero_cnt = -1;
    if (shortcyc[0] == 0 && zero_cnt == -1)
        current = 0;
    else
        current = extremes[2 * shortcyc[0] + 1 + zero_cnt]->i;
    for (i = 0; i < npaths - 1; i++)
    {
        if (shortcyc[i + 1] == 0 && zero_cnt == -1)
        {
            arc     = graph_add_arc(graph, current, 0);
            arc->x  = 1.0;
            current = 0;
        }
        else
        {
            arc     = graph_add_arc(graph, current,
                                extremes[2 * shortcyc[i + 1] + zero_cnt]->i);
            arc->x  = 1.0;
            current = extremes[2 * shortcyc[i + 1] + 1 + zero_cnt]->i;
        }
    }
    if (shortcyc[0] == 0 && zero_cnt == -1)
        arc = graph_add_arc(graph, current, 0);
    else
        arc =
        graph_add_arc(graph, current, extremes[2 * shortcyc[0] + zero_cnt]->i);
    arc->x = 1.0;

    free(shortcyc);
    return 0;
}

static void
connect_search(solver_graph *graph, int n, int comp)
{
    int head = 0;
    int *dstack;
    graph_vertex *v, *other;
    graph_arc *arc;

    dstack            = malloc(graph->nv * sizeof(int));
    graph->v[n]->comp = comp;
    dstack[head++]    = n;

    while (head > 0)
    {
        v = graph->v[dstack[--head]];
        for (arc = v->edge; arc; arc = outnext(arc, v))
        {
            other = otherend(arc, v);
            if (!other->comp)
            {
                other->comp    = comp;
                dstack[head++] = other->i;
            }
        }
    }
    free(dstack);
}

static int
sort_extremes(const void *vv, const void *uu)
{
    graph_vertex *v = *(graph_vertex **)vv, *u = *(graph_vertex **)uu;

    if (v->comp < u->comp)
        return -1;
    if (v->comp > u->comp)
        return +1;

    return 0;
}

int
cp_get_xheur_greedy__(void *cp_, void *bac_env_, void *sol_, int *nadded)
{
    int rval;
    cp_prob *cp               = (cp_prob *)cp_;
    cp_exact_bac_env *bac_env = (cp_exact_bac_env *)bac_env_;
    cp_sol *sol               = (cp_sol *)sol_;
    lp_prob *lp               = bac_env->ip->lp;
    cp_sol *tmpsol            = NULL;
    cp_cut *cuts              = NULL;

    tmpsol = cp_create_sol(cp);
    check_null(tmpsol, " out of memory", CLEANUP);

    rval = cp_get_xheur_greedy(cp, bac_env, tmpsol);
    check_rval(rval, "failed", CLEANUP);

    if (bac_env->ip->lowerboundN < sol->val - SOLVER_ZEROPLUS)
        bac_env->ip->lowerboundN = sol->val;
    if (cp->ip->lowerboundG < sol->val - SOLVER_ZEROPLUS)
    {
        cp->ip->lowerboundG = sol->val;
        cp_copy_sol(tmpsol, sol);
    }

    cuts = cp_conv_cut_sol2connect(tmpsol);
    rval = cp_add_lp_cuts(cp, bac_env, &cuts, nadded);

    if (lp->status == SOLVER_LP_SUCCESS &&
        (bac_env->ip->upperboundN <
         bac_env->ip->lowerboundN + bac_env->ip->param->pruning_tol))
    {
        int edge_added;
        rval = cp_add_badvars(cp, bac_env, &edge_added);
    }

CLEANUP:
    cp_free_sol(&tmpsol);
    return rval;
}

int
cp_get_xheur_greedy(cp_prob *cp, cp_exact_bac_env *bac_env, cp_sol *sol)
{
    int rval;
    int i, k;
    graph_vertex **extremes;
    int nextremes, comp;
    graph_arc *lparc, *new, *e, *next;
    graph_vertex *prev, *other;
    solver_graph *graph;
    lp_prob *lp           = bac_env->ip->lp;
    cp_heur_env *heur_env = bac_env->heur;

    rval = stats_start(bac_env->ip->stats->xheur);

    if (lp->sol->graph->nv != cp->n)
        return 0;

    graph = graph_create();
    graph_add_vertices(graph, cp->n);
    for (i = 0; i < cp->n; i++) graph->v[i]->y = lp->sol->graph->v[i]->y;

    graph_arc **edges = malloc(lp->sol->graph->na * sizeof(graph_arc *));

    for (i = 0; i < lp->sol->graph->na; i++) edges[i] = lp->sol->graph->arcs[i];

    for (i = lp->sol->graph->na; i > 0; i--)
    {
        k = rand() % i;
        SWAP(edges[i - 1], edges[k], e);
    }

    qsort(edges, lp->sol->graph->na, sizeof(graph_arc *), sort_edges);

    comp = 1;
    for (k = 0; k < lp->sol->graph->na; k++)
    {
        lparc = edges[k];

        if (lparc->x < 0.3 && graph->na > 3) // TODO: set as param
            break;

        if (graph->v[lparc->tail->i]->deg < 2 &&
            graph->v[lparc->head->i]->deg < 2)
        {
            prev  = graph->v[lparc->tail->i];
            other = NULL;
            for (e = prev->edge; e && (!other || prev->i != lparc->tail->i);
                 e = next)
            {
                other = otherend(e, prev);
                if (other->i != lparc->head->i)
                    next = outnext(e, other) ? outnext(e, other) : other->edge;
                else
                    next = NULL;
                prev = other;
            }
            if (!other || other->i != lparc->head->i)
            {
                new    = graph_add_arc(graph, lparc->tail->i, lparc->head->i);
                new->x = lparc->x;
            }
        }
    }

    /****************************/
    /* Join the extreme points **/
    /****************************/

    comp = 0;
    for (i = 0; i < graph->nv; i++) graph->v[i]->comp = 0;
    if (graph->v[0]->deg == 0)
        graph->v[0]->comp = ++comp;
    for (i = comp; i < graph->nv; i++)
    {
        if (graph->v[i]->deg && graph->v[i]->comp == 0)
        {
            comp++;
            connect_search(graph, i, comp);
        }
    }

    if (comp > 0)
    {
        extremes = malloc(2 * comp * sizeof(graph_vertex *));

        nextremes = 0;
        if (graph->v[0]->deg == 0)
        {
            graph->v[0]->ind      = nextremes;
            extremes[nextremes++] = graph->v[0];
        }
        for (i = 0; i < graph->nv; i++)
        {
            if (graph->v[i]->deg == 1)
            {
                graph->v[i]->ind      = nextremes;
                extremes[nextremes++] = graph->v[i];
            }
        }

        if (nextremes == 2)
        {
            new    = graph_add_arc(graph, extremes[0]->i, extremes[1]->i);
            new->x = 1.0;
        }
        else if (nextremes > 2)
        {
            qsort(extremes, nextremes, sizeof(graph_vertex *), sort_extremes);

            call_random_tour(cp, graph, nextremes, extremes);
        }
        free(extremes);
    }

    cp_get_sol_from_graph(cp, graph, &sol);

    cp_improve_heur_cycle_length(cp, heur_env, sol);
    heur_env->recover_infeas(cp, heur_env, sol);
    heur_env->local_search(cp, heur_env, sol);

    cp->sol_status = SOLVER_FEAS;

    rval = stats_stop(bac_env->ip->stats->xheur, 1);

    graph_free(&graph);
    free(edges);
    return (0);
}
