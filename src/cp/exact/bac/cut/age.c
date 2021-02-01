#include "cp/cp.h"
#include "ip/exact/bac/bac.h"
#include "cp/exact/bac/bac.h"

static int
get_cut_slack_(cp_prob *cp, cp_exact_bac_env *bac_env, double *cut_slack);

int
cp_age_cuts(cp_prob *cp, cp_exact_bac_env *bac_env, int *ndel)
{
    int *del          = NULL;
    double *cut_slack = NULL;
    int rval;
    int i, j;
    cp_cut *cut;
    int nrows, cut_start;

    cp_cut_repo *cuts = bac_env->cuts;
    lp_prob *lp       = bac_env->ip->lp;

    *ndel = 0;

    if (cuts->count)
    {
        cut_slack = malloc(cuts->count * sizeof(double));
        check_null(cut_slack, "Out of memory ", cleanup);
    }
    rval = get_cut_slack_(cp, bac_env, cut_slack);
    check_rval(rval, "get_cut_slack_ failed", cleanup);

    nrows = lp_get_nrows(lp);
    del   = malloc(nrows * sizeof(int));
    check_null(del, "Out of memory ", cleanup);

    for (i = cuts->count - 1; i >= 0; i--)
    {
        cut = cuts->cuts[i];
        if (cut->branch == 0)
        {
            if (cut_slack[i] < SOLVER_IP_BAC_SLACK_DUST)
            {
                cut->age = 0;
            }
            else
            {
                cut->age++;
                if (cut->age == SOLVER_IP_BAC_MAX_CUT_AGE)
                {
                    cp_del_cut_hash(cuts->hash, cut);
                    if (cut->tcount + cut->hcount)
                    {
                        cp_unregister_cut_repo_cliques(cuts, cut);
                        assert(cut->cliques == NULL);
                    }
                    else
                    {
                        cp_free_cut(&cut);
                        cuts->cuts[i] = NULL;
                    }
                }
            }
        }
    }

    if (cp->data->cap)
        cut_start = cp->data->n + 1;
    else
        cut_start = cp->data->n;

    for (i = 0; i < cut_start; i++) del[i] = 0;

    for (i = 0, j = 0; i < cuts->count; i++)
    {
        cut = cuts->cuts[i];
        if (cut &&
            (cut->cliques == cuts->cliques->cliques || cut->logical ||
             cut->cover_edge || cut->cover_vertex || cut->connect || cut->path))
        {
            if (j < i)
                cuts->cuts[j] = cut;
            j++;
            del[cut_start + i] = 0;
        }
        else
        {
            if (cut && cut->tcount + cut->hcount)
                cp_free_cut(&cut);
            del[cut_start + i] = 1;
        }
    }

    if (j < cuts->count)
    {
        rval = lp_del_rows(lp, del);
        check_rval(rval, "lp_del_rows failed", cleanup);
    }

    *ndel       = cuts->count - j;
    cuts->count = j;
    rval        = 0;

cleanup:
    if (cut_slack)
        free(cut_slack);
    if (del)
        free(del);
    return rval;
}

static int
get_cut_slack_(cp_prob *cp, cp_exact_bac_env *bac_env, double *cut_slack)
{
    int rval;
    double *pi  = NULL;
    lp_prob *lp = bac_env->ip->lp;

    cp_cut_repo *cuts = bac_env->cuts;

    int i;
    int ncount, ncuts, nrows;

    ncount = cp->data->n;
    ncuts  = cuts->count;
    if (cp->data->cap)
        nrows = ncount + 1 + ncuts;
    else
        nrows = ncount + ncuts;

    assert(nrows == lp_get_nrows(lp));

    pi = malloc(nrows * sizeof(double));
    check_null(pi, "Out of memory", cleanup);

    rval = lp_get_slack(lp, pi);
    check_rval(rval, "lp_get_slack failed", cleanup);
    check_null(pi, "lp_get_slack failed", cleanup);

    if (cut_slack)
    {
        for (i = 0; i < ncuts; i++)
        {
            if (cuts->cuts[i]->sense == 'G')
                cut_slack[i] = -pi[ncount + 1 + i];
            else
                cut_slack[i] = pi[ncount + 1 + i];
        }
    }

    rval = 0;
cleanup:
    free(pi);
    return rval;
}

int
cp_age_vars(cp_prob *cp, cp_exact_bac_env *bac_env, int *ndel)
{
    int *del = NULL;
    int rval, ncols;
    graph_arc *arc;
    int i;
    double *x  = NULL;
    double *rc = NULL;

    lp_prob *lp = bac_env->ip->lp;

    solver_graph *graph = lp->graph;

    *ndel = 0;

    ncols = lp_get_ncols(lp);

    x = malloc(ncols * sizeof(double));
    check_null(x, "Out of memory", cleanup);
    rval = lp_get_x(lp, x);
    check_rval(rval, "lp_get_x failed", cleanup);

    rc = malloc(graph->na * sizeof(double));
    check_null(rc, "Out of memory ", cleanup);
    rval = lp_get_rc(lp, rc, graph->nv, graph->nv + graph->na - 1);
    check_rval(rval, "lp_get_rc failed", cleanup);

    *ndel = 0;
    del   = malloc(ncols * sizeof(int));
    check_null(del, "Out of memory ", cleanup);
    for (i = 0; i < ncols; i++) del[i] = 0;

    for (i = graph->na - 1; i >= 0; i--)
    {
        arc = graph->arcs[i];
        if (arc->fixed == 0 && arc->branch == 0)
        {
            if (x[graph->nv + arc->ind] > SOLVER_IP_BAC_VAR_DUST)
            {
                arc->age = 0;
            }
            else
            {
                arc->age++;
                if (arc->age >= SOLVER_IP_BAC_MAX_VAR_AGE)
                {
                    del[graph->nv + arc->ind] = 1;
                    (*ndel)++;
                    graph_del_arc(graph, &arc);
                }
            }
        }
    }

    if (*ndel)
    {
        rval = lp_del_cols(lp, del);
        check_rval(rval, "", cleanup);
    }

    rval = 0;

cleanup:
    if (del)
        free(del);
    if (x)
        free(x);
    if (rc)
        free(rc);
    return rval;
}
