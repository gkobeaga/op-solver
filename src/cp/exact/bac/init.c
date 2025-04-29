#include "ip/ip.h"
#include "ip/exact/bac/bac.h"
#include "cp/cp.h"
#include "cp/heur/heur.h"
#include "cp/exact/exact.h"
#include "cp/exact/bac/bac.h"

static int
load_first_lp(cp_prob *cp, cp_exact_bac_env *bac_env);

int
cp_init_exact_bac(cp_prob *cp, cp_exact_bac_env *bac_env)
{
    int rval = 0;
    int nadded;

    ip_exact_bac_env *ip = bac_env->ip;
    lp_prob *lp          = bac_env->ip->lp;

    ip->orig_env  = bac_env;
    ip->orig_prob = cp;

    solver_data *data = cp->data;
    graph_arc *arc;

    /**************************/
    /* Build and Configure IP */
    /**************************/

    bac_env->ip->infeas_bound = cp->data->tot_obj_node;
    cp->ip->upperboundG       = cp->data->tot_obj_node;

    bac_env->cuts = cp_create_cut_repo(cp);

    if (cp->sol->val > cp->ip->lowerboundG + SOLVER_ZEROPLUS)
    {
        if (bac_env->ip->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("Setting lowerbound to the initial bound: %.2f\n",
                   cp->sol->val);
        }
        cp->ip->lowerboundG = cp->sol->val;
    }
    if (cp->sol->val > ip->lowerboundN + SOLVER_ZEROPLUS)
    {
        ip->lowerboundN = cp->sol->val;
    }

    /*******************************************/
    /* Create the graph of the First LP (k-NN) */
    /*******************************************/
    graph_add_vertices(lp->graph, cp->n);
    for (int i = 0; i < cp->n; i++)
    {
        lp->graph->v[i]->obj = cp->data->obj_node[i];
    }

    data_get_k_nearest(data, 6);

    for (int i = 0; i < data->map->kn_ecount; i++)
    {
        arc       = graph_add_arc(lp->graph, data->map->kn_elist[2 * i],
                                  data->map->kn_elist[2 * i + 1]);
        arc->cost = data_get_norm(data, data->map->kn_elist[2 * i],
                                  data->map->kn_elist[2 * i + 1]);
    }

    rval = load_first_lp(cp, bac_env);
    check_rval(rval, "failed\n", CLEANUP);

    /**************************/
    /* Optimize first LP      */
    /**************************/
    rval = stats_start(bac_env->stats->lp_opt);
    check_rval(rval, "stats_start failed\n", CLEANUP);
    rval = lp_opt_dual(lp);
    check_rval(rval, "lp_opt failed\n", CLEANUP);
    rval = stats_stop(bac_env->stats->lp_opt, 0);
    check_rval(rval, "stats_stop failed\n", CLEANUP);

    if (rval == SOLVER_LP_INFEASIBLE)
    {
        printf("Initial lp infeasible\n");
        goto CLEANUP;
    }

    rval = cp_update_lp_sol(cp, bac_env);
    check_rval(rval, "failed", CLEANUP);

    /**************************/
    /* Add bad vars to LP     */
    /**************************/
    rval = ip->pricing_loop(cp, bac_env, &nadded);
    check_rval(rval, "", CLEANUP);
    cp->ip->upperboundG = ip->upperboundN;

    if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("Initial LP value: %.6f\n", lp->sol->val);
    }

CLEANUP:

    return rval;
}

int
load_first_lp(cp_prob *cp, cp_exact_bac_env *bac_env)
{
    int rval = 0;
    lp_data *lp_data;
    lp_prob *lp = bac_env->ip->lp;

    if (bac_env->ip->verbosity >= SOLVER_VERBOSITY_INFO)
        printf("Loading lp...");

    lp->graph->v[cp->data->from]->fixed  = 1;
    lp->graph->v[cp->data->from]->branch = 0;
    lp->graph->fixed                     = lp->graph->v[cp->data->from];
    if (cp->data->from != cp->data->to)
    {
        lp->graph->v[cp->data->to]->fixed        = 1;
        lp->graph->v[cp->data->to]->branch       = 0;
        lp->graph->v[cp->data->to]->fixed_prev   = lp->graph->v[cp->data->from];
        lp->graph->v[cp->data->from]->fixed_next = lp->graph->v[cp->data->to];
    }

    // Lower and upper bounds of variables
    lp_data =
    cp_build_lp_data(cp, bac_env, 0, lp->graph->nv + lp->graph->na - 1);
    check_null(lp_data, "", CLEANUP);

    // Distance limitation constraint
    lp_data->rhs[cp->n]   = cp->data->cap;
    lp_data->sense[cp->n] = 'L';

    rval = lp_init(lp, lp_data);
    check_rval(rval, " failed", CLEANUP);

    if (bac_env->ip->verbosity >= SOLVER_VERBOSITY_INFO)
    {
        printf("LP has:  %d rows  %d columns  %d nonzeros\n", lp_get_nrows(lp),
               lp_get_ncols(lp), lp_get_nnonzeros(lp));
    }

CLEANUP:
    lp_free_data(&lp_data);

    return rval;
}
