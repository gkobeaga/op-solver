#include "cp/cp.h"
#include "kp/kp.h"

int
cp_sep_cover_vertex(cp_prob *cp, cp_exact_bac_env *bac_env, int *cutcount,
                    cp_cut **cuts)
{
    int rval = 0;
    int i, nfrac, nverts = 0;
    int nint = 0;
    double val;
    double capacity;
    double profit_sum;
    int *verts = NULL;
    cp_cut *cut;
    cp_cut_cover_vertex *cover;
    solver_data *data = NULL;
    kp_prob *kp       = NULL;

    solver_graph *graph = bac_env->ip->lp->sol->graph;

    *cutcount = 0;
    *cuts     = NULL;

    graph_vertex *v;
    int *intind  = NULL;
    int *fracind = NULL;

    verts = malloc(graph->nv * sizeof(int));

    intind = malloc(graph->nv * sizeof(int));

    capacity = -cp->ip->upperboundG - 1.0;
    for (i = 0, nint = 0, nfrac = 0; i < graph->nv; i++)
    {
        v = graph->v[i];
        if (v->y >= SOLVER_ONEMINUS)
        {
            intind[nint++] = i;
            capacity += v->obj;
        }
        else if (v->deg)
        {
            nfrac++;
            capacity += v->obj;
        }
    }

    if (capacity < 0)
        goto CLEANUP;

    if (nfrac == 0)
    {
        printf("The given vector is integer\n");
        assert(capacity != cp->ip->upperboundG);
        goto CLEANUP;
    }

    assert(nfrac == graph->n3v - nint);
    nfrac     = graph->n3v - nint;
    data      = data_create();
    data->cap = capacity;
    data->n   = nfrac;
    kp        = kp_create_prob(data);
    kp->cap   = capacity;

    fracind    = malloc(nfrac * sizeof(int));
    profit_sum = 0.0;
    for (i = 0, nfrac = 0; i < graph->nv; i++)
    {
        v = graph->v[i];
        if (v->deg && v->y < SOLVER_ONEMINUS)
        {
            if (kp->cap > v->obj)
            {
                kp->w[nfrac] = v->obj;
                kp->p[nfrac] = (1.0 - v->y);
                profit_sum += (1.0 - v->y);
                fracind[nfrac++] = i;
                capacity -= v->obj;
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
        for (i = 0; i < data->n; i++)
        {
            val += kp->p[i];
            kp->sol->selected[i] = 1;
        }
        if (val > SOLVER_ZEROPLUS)
        {
            for (i = 0; i < nfrac; i++) verts[nverts++] = fracind[i];
        }
        else
            goto CLEANUP;
    }
    else
    {
        kp_opt(kp, bac_env->kp, kp->sol);

        for (i = 0; i < data->n; i++)
        {
            if (kp->sol->selected[i])
                val += kp->p[i];
        }
        if (val > SOLVER_ZEROPLUS)
        {

            kp_get_mincover(kp, kp->sol);
            for (i = 0, nverts = 0; i < nfrac; i++)
            {
                if (kp->sol->selected[i] == 0)
                {
                    verts[nverts++] = fracind[i];
                }
            }
        }
    }

    if (val > SOLVER_ZEROPLUS)
    {
        cut = cp_create_cut();
        check_null(cut, "out of memory", CLEANUP)

        for (i = 0; i < nint; i++) verts[nverts++] = intind[i];
        assert(nverts == nint + (nfrac - kp->sol->ns));

        cover = malloc(sizeof(cp_cut_cover_vertex));
        clique_conv_array2clique(verts, nverts, &cover->verts);

        cut->sense        = 'L';
        cut->rhs          = (double)nverts - 1.0;
        cut->cover_vertex = cover;

        cut->next = *cuts;
        *cuts     = cut;
        (*cutcount)++;
    }

CLEANUP:
    if (data)
        data_free(&data);
    if (kp)
        kp_free_prob(&kp);
    if (verts)
        free(verts);
    if (intind)
        free(intind);
    if (fracind)
        free(fracind);
    return rval;
}
