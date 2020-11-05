#include "op-solver.h"

#define LOOP_FULL_IN (500) /* to force a full price after 25 inner loops */

static int
sep_loop_inner_(cp_prob *cp, cp_exact_bac_env *bac_env, int *cut_added,
                int *edge_added),
sep_loop_middle_(cp_prob *cp, cp_exact_bac_env *bac_env, int *cut_added,
                 int *edge_added),
sep_loop_middle_2_(cp_prob *cp, cp_exact_bac_env *bac_env, int *cut_added,
                   int *edge_added),
sep_loop_outer_(cp_prob *cp, cp_exact_bac_env *bac_env, int *cut_added,
                int *edge_added),
primal_heur_(cp_prob *cp, cp_exact_bac_env *bac_env, int *cut_added);

int
cp_sep_loop(void *_cp, void *_bac_env, int *cut_added)
{
    int rval = 0;
    int edge_added;

    cp_prob *cp               = (cp_prob *)_cp;
    cp_exact_bac_env *bac_env = (cp_exact_bac_env *)_bac_env;

    lp_prob *lp = bac_env->ip->lp;

    *cut_added                    = 0;
    cp_exact_bac_stats *bac_stats = bac_env->stats;

    rval = stats_start(bac_stats->sep_loop);
    check_rval(rval, "stats_start failed", cleanup);

    if (bac_env->param->sep_loop == CP_SEP_LOOP_BI_LEVEL)
        sep_loop_middle_2_(cp, bac_env, cut_added, &edge_added);
    else if (bac_env->param->sep_loop == CP_SEP_LOOP_THREE_LEVEL)
        sep_loop_outer_(cp, bac_env, cut_added, &edge_added);

cleanup:

    if (!rval && lp->status == SOLVER_LP_INFEASIBLE)
    {
        printf("LP is infeasible in sep_loop\n");
    }

    stats_stop(bac_stats->sep_loop, !rval);

    return rval;
}

#define MIDDLE_LOOP_APPLY_INNER()                                              \
    rval = sep_loop_inner_(cp, bac_env, &cut_added_, &edge_added_);            \
    *cut_added += edge_added_;                                                 \
    *edge_added += cut_added_;                                                 \
    if (lp->status == SOLVER_LP_INFEASIBLE)                                    \
    {                                                                          \
        stats_stop(bac_stats->sep_loop_middle, 1);                             \
        stats_stop_if_active(bac_stats->sep_loop_middle_it, 0);                \
        goto cleanup;                                                          \
    }                                                                          \
    if (bac_env->ip->upperboundN <                                             \
        bac_env->ip->lowerboundN + bac_param->pruning_tol)                     \
    {                                                                          \
        stats_stop(bac_stats->sep_loop_middle, 0);                             \
        stats_stop_if_active(bac_stats->sep_loop_middle_it, 0);                \
        goto cleanup;                                                          \
    }                                                                          \
    if (stats_get_current_time(bac_stats->total) > bac_param->time_limit)      \
    {                                                                          \
        stats_stop(bac_stats->sep_loop_middle, 0);                             \
        stats_stop_if_active(bac_stats->sep_loop_middle_it, 0);                \
        goto cleanup;                                                          \
    }                                                                          \
    if (lp->sol->integral && lp->graph->connected)                             \
    {                                                                          \
        stats_stop(bac_stats->sep_loop_middle, 1);                             \
        stats_stop_if_active(bac_stats->sep_loop_middle_it, 0);                \
        goto cleanup;                                                          \
    }

#define MIDDLE_LOOP_SEP_CUT(cut)                                               \
    if (bac_param->sep_##cut)                                                  \
    {                                                                          \
        rval = stats_start(bac_stats->sep_##cut);                              \
        check_rval(rval, "stats_start failed", cleanup);                       \
        rval = cp_sep_##cut(cp, bac_env, &cutcount, &cuts);                    \
        check_rval(rval, "separation", cleanup);                               \
                                                                               \
        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)                       \
            printf("Found %2d %s cuts in %.2ld seconds\n", cutcount,           \
                   bac_stats->sep_##cut->name,                                 \
                   stats_get_current_time(bac_stats->sep_##cut));              \
                                                                               \
        if (cutcount)                                                          \
        {                                                                      \
            rval = cp_add_lp_cuts(cp, bac_env, &cuts, &cut_added_);            \
            check_rval(rval, "failed", cleanup);                               \
            stats_stop(bac_stats->sep_##cut, cut_added_);                      \
                                                                               \
            if (lp->status == SOLVER_LP_INFEASIBLE)                            \
            {                                                                  \
                stats_stop(bac_stats->sep_loop_middle, 0);                     \
                stats_stop(bac_stats->sep_loop_middle_it, 0);                  \
                goto cleanup;                                                  \
            }                                                                  \
                                                                               \
            *cut_added += cut_added_;                                          \
                                                                               \
            primal_heur_(cp, bac_env, &cut_added_);                            \
            *cut_added += cut_added_;                                          \
                                                                               \
            if (lp->status == SOLVER_LP_INFEASIBLE ||                          \
                cp->ip->lowerboundG >= cp->ip->upperboundG)                    \
            {                                                                  \
                stats_stop(bac_stats->sep_loop_middle, 0);                     \
                stats_stop(bac_stats->sep_loop_middle_it, 0);                  \
                goto cleanup;                                                  \
            }                                                                  \
            rval = cp_add_badvars(cp, bac_env, &edge_added_);                  \
            check_rval(rval, "failed", cleanup);                               \
            *edge_added += edge_added_;                                        \
            edge_added_it += edge_added_;                                      \
            if (bac_env->ip->upperboundN <                                     \
                bac_env->ip->lowerboundN + bac_param->pruning_tol)             \
            {                                                                  \
                stats_stop(bac_stats->sep_loop_middle, 1);                     \
                stats_stop(bac_stats->sep_loop_middle_it, 0);                  \
                goto cleanup;                                                  \
            }                                                                  \
            if (lp->sol->integral && lp->sol->graph->connected)                \
            {                                                                  \
                rval = cp_add_badvars(cp, bac_env, &edge_added_);              \
                check_rval(rval, "failed", cleanup);                           \
                *edge_added += edge_added_;                                    \
                edge_added_it += edge_added_;                                  \
                if (!edge_added_)                                              \
                {                                                              \
                    stats_stop(bac_stats->sep_loop_middle, 1);                 \
                    stats_stop(bac_stats->sep_loop_middle_it, 0);              \
                    goto cleanup;                                              \
                }                                                              \
            }                                                                  \
            MIDDLE_LOOP_APPLY_INNER();                                         \
        }                                                                      \
        else                                                                   \
            rval = stats_stop(bac_stats->sep_##cut, 0);                        \
    }

static int
sep_loop_middle_(cp_prob *cp, cp_exact_bac_env *bac_env, int *cut_added,
                 int *edge_added)
{
    int rval = 0;
    int cutcount, cut_added_, edge_added_, cut_added_it, edge_added_it;
    double oldval;
    cp_cut *cuts  = NULL;
    int loopcount = 0;
    lp_prob *lp   = bac_env->ip->lp;

    cp_exact_bac_stats *bac_stats = bac_env->stats;
    cp_exact_bac_param *bac_param = bac_env->param;

    int oldsep     = bac_env->param->sep_sec_exact;
    int oldsrk     = bac_env->param->srk_rule;
    int olds2      = bac_env->param->srk_s2;
    int olds3      = bac_env->param->srk_s3;
    int oldmaxv    = bac_env->param->sec_max_vin;
    int oldmaxviol = bac_env->param->sec_max_viol;

    *cut_added  = 0;
    *edge_added = 0;

    rval = stats_start(bac_stats->sep_loop_middle);
    check_rval(rval, "stats_start failed", cleanup);

    edge_added_it = 0;
    MIDDLE_LOOP_APPLY_INNER();

    do
    {
        rval = stats_start(bac_stats->sep_loop_middle_it);
        check_rval(rval, "failed", cleanup);

        oldval        = lp->sol->val;
        cut_added_it  = 0;
        edge_added_it = 0;

        MIDDLE_LOOP_SEP_CUT(blossom_mst);

        MIDDLE_LOOP_SEP_CUT(blossom_fast);

        MIDDLE_LOOP_SEP_CUT(blossom_ghfast);

        MIDDLE_LOOP_SEP_CUT(sec_exact);

        MIDDLE_LOOP_SEP_CUT(cover_cycle);

        if (bac_param->sep_sec_cc_2)
        {
            bac_env->param->sep_sec_exact = CP_SEC_HONG_TS;
            bac_env->param->srk_s3        = 0;
            bac_env->param->sec_max_vin = bac_env->param->sec_max_vout = 1;
            MIDDLE_LOOP_SEP_CUT(sec_exact);
            bac_env->param->sep_sec_exact = oldsep;
            bac_env->param->srk_rule      = oldsrk;
            bac_env->param->srk_s2        = olds2;
            bac_env->param->srk_s3        = olds3;
            bac_env->param->sec_max_vout  = oldmaxv;
            bac_env->param->sec_max_vin   = oldmaxv;
            bac_env->param->sec_max_viol  = oldmaxviol;
        }

        rval = stats_stop(bac_stats->sep_loop_middle_it, 0);
        check_rval(rval, "stats_stop failed", cleanup);

        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("  LOOP vals: newval %f oldval %f\t edge_added %d \t "
                   "loopcount %d \t priceval %f penal %f lowerboundN %.2f "
                   "lowerboundG %.2f\n",
                   lp->sol->val, oldval, edge_added_it, loopcount,
                   bac_env->ip->upperboundN,
                   bac_env->ip->upperboundN - lp->sol->val,
                   bac_env->ip->lowerboundN, cp->ip->lowerboundG);
        }

        loopcount++;

        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("  LOOP MIDDLE %2d: %f  (%.2ld seconds)\n", loopcount,
                   lp->sol->val,
                   stats_get_current_time(bac_stats->sep_loop_middle));
        }

    } while ((edge_added_it || (oldval - lp->sol->val) / oldval * 100 > 1) &&
             loopcount < LOOP_FULL_IN);

    rval = stats_stop(bac_stats->sep_loop_middle, 0);
    check_rval(rval, "stats_stop failed", cleanup);

cleanup:

    bac_env->param->sep_sec_exact = oldsep;
    bac_env->param->srk_rule      = oldsrk;
    bac_env->param->srk_s2        = olds2;
    bac_env->param->srk_s3        = olds3;
    bac_env->param->sec_max_vout  = oldmaxv;
    bac_env->param->sec_max_vin   = oldmaxv;
    bac_env->param->sec_max_viol  = oldmaxviol;

    if (!rval && lp->status == SOLVER_LP_INFEASIBLE)
    {
        printf("LP is infeasible in middle sep_loop\n");
    }

    return rval;
}

#define OUTER_LOOP_APPLY_MIDDLE()                                              \
    rval = sep_loop_middle_(cp, bac_env, &cut_added_, &edge_added_);           \
    *edge_added += edge_added_;                                                \
    *cut_added += cut_added_;                                                  \
    if (lp->status == SOLVER_LP_INFEASIBLE)                                    \
    {                                                                          \
        stats_stop(bac_stats->sep_loop_outer, 1);                              \
        stats_stop_if_active(bac_stats->sep_loop_outer_it, 0);                 \
        goto cleanup;                                                          \
    }                                                                          \
    if (bac_env->ip->upperboundN <                                             \
        bac_env->ip->lowerboundN + bac_param->pruning_tol)                     \
    {                                                                          \
        stats_stop(bac_stats->sep_loop_outer, 0);                              \
        stats_stop_if_active(bac_stats->sep_loop_outer_it, 0);                 \
        goto cleanup;                                                          \
    }                                                                          \
    if (stats_get_current_time(bac_stats->total) > bac_param->time_limit)      \
    {                                                                          \
        stats_stop(bac_stats->sep_loop_outer, 0);                              \
        stats_stop_if_active(bac_stats->sep_loop_outer_it, 0);                 \
        goto cleanup;                                                          \
    }                                                                          \
    if (lp->sol->integral && lp->graph->connected)                             \
    {                                                                          \
        stats_stop(bac_stats->sep_loop_outer, 1);                              \
        stats_stop_if_active(bac_stats->sep_loop_outer_it, 0);                 \
        goto cleanup;                                                          \
    }

#define OUTER_LOOP_SEP_CUT(cut)                                                \
    if (bac_param->sep_##cut)                                                  \
    {                                                                          \
                                                                               \
        rval = stats_start(bac_stats->sep_##cut);                              \
        check_rval(rval, "stats_start failed", cleanup);                       \
        rval = cp_sep_##cut(cp, bac_env, &cutcount, &cuts);                    \
        check_rval(rval, "separation failed", cleanup);                        \
                                                                               \
        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)                       \
            printf("Found %2d %s cuts in %.2ld seconds\n", cutcount,           \
                   bac_stats->sep_##cut->name,                                 \
                   stats_get_current_time(bac_stats->sep_##cut));              \
                                                                               \
        if (cutcount)                                                          \
        {                                                                      \
            rval = cp_add_lp_cuts(cp, bac_env, &cuts, &cut_added_);            \
            check_rval(rval, "failed", cleanup);                               \
            rval = stats_stop(bac_stats->sep_##cut, cut_added_);               \
            *cut_added += cut_added_;                                          \
                                                                               \
            if (lp->status == SOLVER_LP_INFEASIBLE)                            \
            {                                                                  \
                stats_stop(bac_stats->sep_loop_outer, 1);                      \
                stats_stop(bac_stats->sep_loop_outer_it, 0);                   \
                goto cleanup;                                                  \
            }                                                                  \
                                                                               \
            if (lp->sol->integral && lp->sol->graph->connected)                \
            {                                                                  \
                rval = cp_add_badvars(cp, bac_env, &edge_added_);              \
                *edge_added += edge_added_;                                    \
                edge_added_it += edge_added_;                                  \
                check_rval(rval, "failed", cleanup);                           \
                if (!edge_added_)                                              \
                {                                                              \
                    stats_stop(bac_stats->sep_loop_outer, 1);                  \
                    stats_stop(bac_stats->sep_loop_outer_it, 0);               \
                    goto cleanup;                                              \
                }                                                              \
            }                                                                  \
            primal_heur_(cp, bac_env, &cut_added_);                            \
            *cut_added += cut_added_;                                          \
            if (lp->status == SOLVER_LP_INFEASIBLE ||                          \
                cp->ip->lowerboundG >= cp->ip->upperboundG)                    \
            {                                                                  \
                stats_stop(bac_stats->sep_loop_outer, 1);                      \
                stats_stop(bac_stats->sep_loop_outer_it, 0);                   \
                goto cleanup;                                                  \
            }                                                                  \
            rval = cp_add_badvars(cp, bac_env, &edge_added_);                  \
            check_rval(rval, "failed", cleanup);                               \
            *edge_added += edge_added_;                                        \
            edge_added_it += edge_added_;                                      \
            if (bac_env->ip->upperboundN <                                     \
                bac_env->ip->lowerboundN + bac_param->pruning_tol)             \
            {                                                                  \
                stats_stop(bac_stats->sep_loop_outer, 1);                      \
                stats_stop(bac_stats->sep_loop_outer_it, 0);                   \
                goto cleanup;                                                  \
            }                                                                  \
            if (lp->sol->integral && lp->sol->graph->connected)                \
            {                                                                  \
                rval = cp_add_badvars(cp, bac_env, &edge_added_);              \
                *edge_added += edge_added_;                                    \
                edge_added_it += edge_added_;                                  \
                check_rval(rval, "failed", cleanup);                           \
                if (!edge_added_)                                              \
                {                                                              \
                    stats_stop(bac_stats->sep_loop_outer, 1);                  \
                    stats_stop(bac_stats->sep_loop_outer_it, 0);               \
                    goto cleanup;                                              \
                }                                                              \
            }                                                                  \
            OUTER_LOOP_APPLY_MIDDLE();                                         \
        }                                                                      \
        else                                                                   \
            rval = stats_stop(bac_stats->sep_##cut, 0);                        \
    }

#define LOOP_FULL_OUT (200)
static int
sep_loop_outer_(cp_prob *cp, cp_exact_bac_env *bac_env, int *cut_added,
                int *edge_added)
{
    int rval = 0;
    int cutcount, cut_added_, edge_added_, edge_added_it;
    double oldval;
    cp_cut *cuts   = NULL;
    int loopcount  = 0;
    cp_sol *tmpsol = NULL;
    lp_prob *lp    = bac_env->ip->lp;

    *edge_added                   = 0;
    cp_exact_bac_stats *bac_stats = bac_env->stats;
    cp_exact_bac_param *bac_param = bac_env->param;

    rval = stats_start(bac_stats->sep_loop_outer);
    check_rval(rval, "stats_start failed", cleanup);

    OUTER_LOOP_APPLY_MIDDLE();

    loopcount = 0;
    do
    {
        oldval        = lp->sol->val;
        edge_added_it = 0;

        rval = stats_start(bac_stats->sep_loop_outer_it);
        check_rval(rval, "stats_start failed", cleanup);

        oldval = lp->sol->val;

        OUTER_LOOP_SEP_CUT(cover_edge);

        OUTER_LOOP_SEP_CUT(cover_vertex);

        OUTER_LOOP_SEP_CUT(path);

        loopcount++;

        rval = stats_stop(bac_stats->sep_loop_outer_it, 1);
        check_rval(rval, "stats_stop failed", cleanup);

        if (stats_get_current_time(bac_stats->total) > bac_param->time_limit)
            goto cleanup;

        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("  LOOP OUT %2d: %f (impr %.3f) (%.2ld seconds)\n",
                   loopcount, lp->sol->val,
                   (oldval - lp->sol->val) / oldval * 100,
                   stats_get_current_time(bac_stats->sep_loop));
        }

    } while (
    ((edge_added_it) || ((oldval - lp->sol->val) / oldval * 100 > 1)) &&
    loopcount < LOOP_FULL_OUT &&
    bac_env->ip->upperboundN >
    bac_env->ip->lowerboundN + bac_param->pruning_tol);

    rval = stats_stop(bac_stats->sep_loop_outer, 0);
    check_rval(rval, "stats_stop failed", cleanup);

cleanup:

    if (!rval && lp->status == SOLVER_LP_INFEASIBLE)
    {
        printf("LP is infeasible in outer sep_loop\n");
    }

    return rval;
}

static int
sep_loop_middle_2_(cp_prob *cp, cp_exact_bac_env *bac_env, int *cut_added,
                   int *edge_added)
{
    int rval = 0;
    int cutcount, cut_added_, edge_added_, cut_added_it, edge_added_it;
    double oldval;
    cp_cut *cuts  = NULL;
    int loopcount = 0;
    lp_prob *lp   = bac_env->ip->lp;

    cp_exact_bac_stats *bac_stats = bac_env->stats;
    cp_exact_bac_param *bac_param = bac_env->param;

    int oldsep     = bac_env->param->sep_sec_exact;
    int oldsrk     = bac_env->param->srk_rule;
    int olds2      = bac_env->param->srk_s2;
    int olds3      = bac_env->param->srk_s3;
    int oldmaxv    = bac_env->param->sec_max_vin;
    int oldmaxviol = bac_env->param->sec_max_viol;

    *cut_added  = 0;
    *edge_added = 0;

    rval = stats_start(bac_stats->sep_loop_middle);
    check_rval(rval, "stats_start failed", cleanup);

    edge_added_it = 0;
    MIDDLE_LOOP_APPLY_INNER();

    do
    {
        rval = stats_start(bac_stats->sep_loop_middle_it);
        check_rval(rval, "failed", cleanup);

        oldval        = lp->sol->val;
        cut_added_it  = 0;
        edge_added_it = 0;

        MIDDLE_LOOP_SEP_CUT(blossom_mst);

        MIDDLE_LOOP_SEP_CUT(blossom_fast);

        MIDDLE_LOOP_SEP_CUT(blossom_ghfast);

        MIDDLE_LOOP_SEP_CUT(sec_exact);

        MIDDLE_LOOP_SEP_CUT(cover_cycle);

        MIDDLE_LOOP_SEP_CUT(cover_edge);

        MIDDLE_LOOP_SEP_CUT(cover_vertex);

        MIDDLE_LOOP_SEP_CUT(path);

        rval = stats_stop(bac_stats->sep_loop_middle_it, 0);
        check_rval(rval, "stats_stop failed", cleanup);

        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("  LOOP vals: newval %f oldval %f\t edge_added %d \t "
                   "loopcount %d \t priceval %f penal %f lowerboundN %.2f "
                   "lowerboundG %.2f\n",
                   lp->sol->val, oldval, edge_added_it, loopcount,
                   bac_env->ip->upperboundN,
                   bac_env->ip->upperboundN - lp->sol->val,
                   bac_env->ip->lowerboundN, cp->ip->lowerboundG);
        }

        loopcount++;

        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
        {
            printf("  LOOP MIDDLE %2d: %f  (%.2ld seconds)\n", loopcount,
                   lp->sol->val,
                   stats_get_current_time(bac_stats->sep_loop_middle));
        }

    } while ((edge_added_it || (oldval - lp->sol->val) / oldval * 100 > 1) &&
             loopcount < LOOP_FULL_IN);

    rval = stats_stop(bac_stats->sep_loop_middle, 0);
    check_rval(rval, "stats_stop failed", cleanup);

cleanup:

    bac_env->param->sep_sec_exact = oldsep;
    bac_env->param->srk_rule      = oldsrk;
    bac_env->param->srk_s2        = olds2;
    bac_env->param->srk_s3        = olds3;
    bac_env->param->sec_max_vout  = oldmaxv;
    bac_env->param->sec_max_vin   = oldmaxv;
    bac_env->param->sec_max_viol  = oldmaxviol;

    if (!rval && lp->status == SOLVER_LP_INFEASIBLE)
    {
        printf("LP is infeasible in middle sep_loop\n");
    }

    return rval;
}

static int
sep_loop_inner_(cp_prob *cp, cp_exact_bac_env *bac_env, int *cut_added,
                int *edge_added)
{
    int rval = 0;
    int cutcount, cut_added_, cut_added_it;
    int edge_added_ = 0;
    double newval = 0.0, oldval = 0.0;
    cp_cut *cuts;
    cp_sol *sol    = cp->sol;
    cp_sol *tmpsol = NULL;
    lp_prob *lp    = bac_env->ip->lp;

    cp_exact_bac_stats *bac_stats = bac_env->stats;
    cp_exact_bac_param *bac_param = bac_env->param;

    *edge_added = 0;

    rval = stats_start(bac_stats->sep_loop_inner);
    check_rval(rval, "stats_start failed", cleanup);

restart_loop:

    do
    {
        do
        {
            rval = stats_start(bac_stats->sep_loop_inner_it);
            check_rval(rval, "stats_start failed", cleanup);

            oldval = lp->sol->val;
            newval = oldval;

            cut_added_   = 0;
            cut_added_it = 0;

            rval = stats_start(bac_stats->sep_sec_comps);
            check_rval(rval, "stats_start failed", cleanup);
            rval = cp_sep_sec_comps(cp, bac_env, &cutcount, &cuts);
            check_rval(rval, "failed", cleanup);

            if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
                printf("Found %2d %s cuts in %.2ld seconds\n", cutcount,
                       bac_stats->sep_sec_comps->name,
                       stats_get_current_time(bac_stats->sep_sec_comps));

            if (cutcount)
            {
                rval = cp_add_lp_cuts(cp, bac_env, &cuts, &cut_added_);
                check_rval(rval, "failed", cleanup);
                stats_stop(bac_stats->sep_sec_comps, cut_added_);

                *cut_added += cut_added_;
                cut_added_it += cut_added_;

                if (lp->status == SOLVER_LP_INFEASIBLE)
                {
                    rval = stats_stop(bac_stats->sep_loop_inner_it, 1);
                    goto cleanup;
                }
            }
            else
                stats_stop(bac_stats->sep_sec_comps, 0);

            if (bac_param->sep_logical)
            {
                rval = stats_start(bac_stats->sep_logical);
                check_rval(rval, "stats_start failed", cleanup);
                rval = cp_sep_logical(cp, bac_env, &cutcount, &cuts);
                check_rval(rval, "", cleanup);

                if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
                    printf("Found %2d %s cuts in %.2ld seconds\n", cutcount,
                           bac_stats->sep_logical->name,
                           stats_get_current_time(bac_stats->sep_logical));

                if (cutcount)
                {
                    rval = cp_add_lp_cuts(cp, bac_env, &cuts, &cut_added_);
                    check_rval(rval, "failed", cleanup);
                    rval = stats_stop(bac_stats->sep_logical, cut_added_);

                    *cut_added += cut_added_;
                    cut_added_it += cut_added_;

                    if (lp->status == SOLVER_LP_INFEASIBLE)
                    {
                        rval = stats_stop(bac_stats->sep_loop_inner_it, 1);
                        goto cleanup;
                    }
                }
                else
                    rval = stats_stop(bac_stats->sep_logical, 0);
            }

            if (stats_get_current_time(bac_stats->total) >
                bac_param->time_limit)
            {
                rval = stats_stop(bac_stats->sep_loop_inner_it, 0);
                goto cleanup;
            }

            edge_added_ = 0;
            if (lp->sol->val <=
                bac_env->ip->lowerboundN + bac_param->pruning_tol)
            {
                rval = cp_add_badvars(cp, bac_env, &edge_added_);
                *edge_added += edge_added_;
                check_rval(rval, "failed", cleanup);
            }

            rval = stats_stop(bac_stats->sep_loop_inner_it,
                              cut_added_it + edge_added_);

            newval = lp->sol->val;

        } while ((newval < oldval - SOLVER_ZEROPLUS || (edge_added_)) &&
                 bac_env->ip->upperboundN >=
                 bac_env->ip->lowerboundN + bac_env->ip->param->pruning_tol);

        cut_added_  = 0;
        edge_added_ = 0;

        tmpsol = cp_create_sol(cp);
        check_null(tmpsol, " out of memory", cleanup);

        rval = stats_start(bac_stats->xheur_sep);
        check_rval(rval, "stats_start failed", cleanup);
        rval = cp_get_xheur_greedy(cp, bac_env, tmpsol);
        check_rval(rval, "failed", cleanup);

        if (tmpsol && tmpsol->val > bac_env->ip->lowerboundN)
        {
            bac_env->ip->lowerboundN = tmpsol->val;
            if (tmpsol->val > cp->ip->lowerboundG)
                cp->ip->lowerboundG = tmpsol->val;
            cp_copy_sol(tmpsol, sol);
            cp->sol_status = SOLVER_FEAS;

            rval = stats_stop(bac_stats->xheur_sep, 1);
        }
        else
        {
            rval = stats_stop(bac_stats->xheur_sep, 0);
        }

        cuts = cp_conv_cut_sol2connect(tmpsol);
        check_null(cuts, "failed", cleanup);
        rval = cp_add_lp_cuts(cp, bac_env, &cuts, &cut_added_);
        check_rval(rval, "failed", cleanup);
        *cut_added += cut_added_;

        cp_free_sol(&tmpsol);

        if (lp->status == SOLVER_LP_INFEASIBLE ||
            cp->ip->lowerboundG >= cp->ip->upperboundG)
            goto cleanup;

        if (lp->sol->val <= bac_env->ip->lowerboundN + bac_param->pruning_tol)
        {
            rval = cp_add_badvars(cp, bac_env, &edge_added_);
            *edge_added += edge_added_;
            check_rval(rval, "failed", cleanup);
        }

    } while ((edge_added_ || cut_added_) &&
             bac_env->ip->upperboundN >=
             bac_env->ip->lowerboundN + bac_env->ip->param->pruning_tol);

    if ((lp->sol->integral && lp->sol->graph->connected))
    {

        if (lp->sol->val > bac_env->ip->lowerboundN + SOLVER_ZEROPLUS)
        {
            cp_get_sol_from_graph(cp, lp->sol->graph, &tmpsol);

            bac_env->ip->lowerboundN = lp->sol->val;
            if (lp->sol->val > cp->ip->lowerboundG + SOLVER_ZEROPLUS)
                cp->ip->lowerboundG = lp->sol->val;
            cp_copy_sol(tmpsol, cp->sol);
            cp_free_sol(&tmpsol);
        }

        rval = cp_add_badvars(cp, bac_env, &edge_added_);
        *edge_added += edge_added_;
        check_rval(rval, "failed", cleanup);
        if (edge_added_)
            goto restart_loop;
    }

cleanup:
    rval = stats_stop(bac_stats->sep_loop_inner, edge_added_);
    return rval;
}

static int
primal_heur_(cp_prob *cp, cp_exact_bac_env *bac_env, int *cut_added)
{
    int rval       = 0;
    cp_sol *tmpsol = NULL;
    cp_cut *cuts   = NULL;

    tmpsol = cp_create_sol(cp);
    check_null(tmpsol, " out of memory", cleanup);

    cp_exact_bac_stats *bac_stats = bac_env->stats;
    lp_prob *lp                   = bac_env->ip->lp;

    *cut_added = 0;

    rval = stats_start(bac_stats->xheur_sep);
    check_rval(rval, "stats_start failed", cleanup);
    rval = cp_get_xheur_greedy(cp, bac_env, tmpsol);
    check_rval(rval, "failed", cleanup);

    if (tmpsol->val > bac_env->ip->lowerboundN + SOLVER_ZEROPLUS)
    {
        rval = stats_stop(bac_stats->xheur_sep, 1);
        check_rval(rval, "failed", cleanup);
        if (bac_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("New lowerbound from primal heuristic: %.2f\n", tmpsol->val);
        bac_env->ip->lowerboundN = tmpsol->val;
        if (tmpsol->val > cp->ip->lowerboundG + SOLVER_ZEROPLUS)
            cp->ip->lowerboundG = tmpsol->val;
        cp_copy_sol(tmpsol, cp->sol);
        cp->sol_status = SOLVER_FEAS;
    }
    else
    {
        rval = stats_stop(bac_stats->xheur_sep, 0);
        check_rval(rval, "failed", cleanup);
    }

    if (cp->ip->lowerboundG >= cp->ip->upperboundG)
        goto cleanup;

    cuts = cp_conv_cut_sol2connect(tmpsol);
    rval = cp_add_lp_cuts(cp, bac_env, &cuts, cut_added);
    if (lp->status == SOLVER_LP_INFEASIBLE)
        goto cleanup;
    check_rval(rval, "failed", cleanup);

cleanup:
    cp_free_sol(&tmpsol);
    return rval;
}
