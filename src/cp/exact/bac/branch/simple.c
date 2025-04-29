#include "op-solver.h"
#include "ip/exact/bac/bac.h"
#include "cp/exact/bac/bac.h"

int
cp_find_branch(void *prob, void *env, ip_branch **bobj)
{
    int rval        = 0;
    graph_arc *edge = NULL;
    *bobj           = NULL;

    cp_prob *cp               = (cp_prob *)prob;
    cp_exact_bac_env *bac_env = (cp_exact_bac_env *)env;

    bac_env->ngot = 0;

    rval = cp_find_branch_edge(cp, bac_env, &edge);
    check_rval(rval, "ip_find_branch failed", CLEANUP);
    check_null(edge, "ip_find_branch failed", CLEANUP);

    *bobj = ip_create_branch();
    check_null(*bobj, "out of memory in ip_find_branch", CLEANUP);

    graph_arc *tedge = graph_find_arc_hash(bac_env->ip->lp->sol->graph->archash,
                                           edge->tail->i, edge->head->i);
    (*bobj)->edge    = edge;

    bac_env->ngot = 1;
    bac_env->ip->branch_count++;

CLEANUP:

    return rval;
}

int
cp_find_branch_edge(cp_prob *cp, cp_exact_bac_env *bac_env, graph_arc **edge)
{
    int rval = 0;
    double maxdiff;
    int i;
    graph_arc *tedge, *best;

    lp_prob *lp = bac_env->ip->lp;

    *edge = NULL;

    best    = NULL;
    maxdiff = -1.0;

    for (i = 0; i < lp->sol->graph->na; i++)
    {
        tedge = lp->sol->graph->arcs[i];
        if (!tedge->fixed && !tedge->branch)
        {
            if (tedge->x < 0.5)
            {
                if (tedge->x > maxdiff)
                {
                    maxdiff = tedge->x;
                    best    = graph_find_arc_hash(lp->graph->archash,
                                                  tedge->tail->i, tedge->head->i);
                }
            }
            else if (tedge->x <= SOLVER_ONEMINUS)
            {

                if ((1.0 - tedge->x) > maxdiff)
                {
                    maxdiff = 1.0 - tedge->x;
                    best    = graph_find_arc_hash(lp->graph->archash,
                                                  tedge->tail->i, tedge->head->i);
                }
            }
        }
    }

    check_null(best, "All edges are either branched or fixed", CLEANUP);

    *edge = best;

CLEANUP:

    return rval;
}
