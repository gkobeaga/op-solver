#include "cp/cp.h"
#include "kp/kp.h"

int
cp_sep_cover_edge(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                  cp_cut **cuts)
{
    int rval = 0;
    int i, k, nfrac;
    int nint = 0;
    double val;
    double capacity;
    double profit_sum;
    cp_cut *cut;
    cp_cut_cover_edge *cover;
    solver_data *data = NULL;
    kp_prob *kp       = NULL;

    solver_graph *graph = bac_env->ip->lp->sol->graph;

    *cutcount = 0;
    *cuts     = NULL;

    graph_arc *arc;
    int *intind  = NULL;
    int *fracind = NULL;

    intind = malloc(graph->na * sizeof(int));

    capacity = -cp->data->cap - 1.0;
    for (i = 0, nint = 0; i < graph->na; i++)
    {
        arc = graph->arcs[i];
        if (arc->x >= SOLVER_ONEMINUS)
            intind[nint++] = i;
        capacity += arc->cost;
    }

    if (graph->na == nint)
    {
        printf("The given vector is integer\n");
        goto CLEANUP;
    }

    nfrac   = graph->na - nint;
    data    = data_create();
    data->n = nfrac;
    kp      = kp_create_prob(data);
    kp->cap = capacity;

    fracind    = malloc(nfrac * sizeof(int));
    profit_sum = 0.0;
    for (i = 0, nfrac = 0; i < graph->na; i++)
    {
        arc = graph->arcs[i];
        if (arc->x < SOLVER_ONEMINUS)
        {
            if (kp->cap > arc->cost && arc->cost >= 1)
            {
                kp->w[nfrac] = arc->cost;
                kp->p[nfrac] = (1.0 - arc->x);
                profit_sum += (1.0 - arc->x);
                fracind[nfrac++] = i;
                capacity -= arc->cost;
            }
            else
            {
                data->n--;
                kp->n--;
                intind[nint++] = i;
            }
        }
    }

    val = 1.0 - profit_sum;
    if (capacity >= 0)
    {
        for (i = 0; i < nfrac; i++)
        {
            val += kp->p[i];
            kp->sol->selected[i] = 1;
        }
        kp->sol->ns = nfrac;
        if (val <= SOLVER_ZEROPLUS)
            goto CLEANUP;
    }
    else
    {
        kp_opt(kp, bac_env->kp, kp->sol);

        for (i = 0; i < nfrac; i++)
        {
            if (kp->sol->selected[i])
                val += kp->p[i];
        }
        if (val > SOLVER_ZEROPLUS)
            kp_get_mincover(kp, kp->sol);
    }

    if (val > SOLVER_ZEROPLUS)
    {
        cut = cp_create_cut();
        check_null(cut, "out of memory", CLEANUP)

        cover = malloc(sizeof(cp_cut_cover_edge));
        cover->arcs =
        malloc(2 * (nint + graph->na - kp->sol->ns) * sizeof(int));

        for (k = 0; k < nint; k++)
        {
            arc                    = graph->arcs[intind[k]];
            cover->arcs[2 * k]     = arc->tail->i;
            cover->arcs[2 * k + 1] = arc->head->i;
        }
        assert(k == nint);

        for (i = 0; i < nfrac; i++)
        {
            if (kp->sol->selected[i] == 0)
            {
                arc                    = graph->arcs[fracind[i]];
                cover->arcs[2 * k]     = arc->tail->i;
                cover->arcs[2 * k + 1] = arc->head->i;
                k++;
            }
        }
        assert(k == nint + (nfrac - kp->sol->ns));

        cover->na       = k;
        cut->sense      = 'L';
        cut->rhs        = (double)k - 1.0;
        cover->strong   = 0;
        cut->cover_edge = cover;

        cut->next = *cuts;
        *cuts     = cut;
        (*cutcount)++;
    }

CLEANUP:
    if (data)
        data_free(&data);
    if (kp)
        kp_free_prob(&kp);
    if (intind)
        free(intind);
    if (fracind)
        free(fracind);
    return rval;
}
