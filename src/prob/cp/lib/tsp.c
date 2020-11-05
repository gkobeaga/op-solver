#include "op-solver.h"

void
cp_conv_sol_cp2tsp(solver_data *data, cp_sol *cpsol, tsp_sol *tspsol)
{
    int i;
    tsp_erase_sol(tspsol);
    for (i = 0; i < cpsol->ns; i++)
    {
        tspsol->cycle[i]     = data->map->fun[cpsol->cycle[i]];
        tspsol->cod_fr[i]    = data->map->fun[cpsol->cod_fr[data->map->inv[i]]];
        tspsol->cod_bk[i]    = data->map->fun[cpsol->cod_bk[data->map->inv[i]]];
        tspsol->selected[i]  = 1;
        tspsol->sposition[i] = i;
    }

    tspsol->tot_n = cpsol->ns;
    tspsol->ns    = cpsol->ns;
    tspsol->val   = cpsol->cap;

    return;
}

void
cp_conv_sol_tsp2cp(solver_data *data, tsp_sol *tspsol, cp_sol *cpsol)
{
    int i;
    data_map *map = data->map;
    cp_erase_sol(cpsol);
    for (i = 0; i < tspsol->tot_n; i++)
    {
        cpsol->cycle[i] = map->inv[tspsol->cycle[i]];

        cpsol->cod_fr[map->inv[tspsol->cycle[i]]] =
        map->inv[tspsol->cod_fr[tspsol->cycle[i]]];

        cpsol->cod_bk[map->inv[tspsol->cycle[i]]] =
        map->inv[tspsol->cod_bk[tspsol->cycle[i]]];

        cpsol->selected[map->inv[tspsol->cycle[i]]] = 1;
        cpsol->val += data->obj_node[map->orig[tspsol->cycle[i]]];
    }
    cpsol->cap   = tspsol->val;
    cpsol->tot_n = map->dom_n;
    cpsol->ns    = map->img_n;
}

int
cp_improve_heur_cycle_length(cp_prob *cp, cp_heur_env *env, cp_sol *sol)
{
    int rval          = 0;
    tsp_prob *tsp     = NULL;
    solver_data *data = cp->data;
    data_map *map     = NULL;

    map = data_emb_map(data, sol->selected);
    tsp = tsp_create_prob(data);

    cp_conv_sol_cp2tsp(data, sol, tsp->sol);

    data_get_k_nearest(data, 10);
    tsp_opt_heur(tsp, env->tsp, tsp->sol);

    cp_conv_sol_tsp2cp(data, tsp->sol, sol);
    data->map = map->prev;

    data_free_map(&map);
    tsp_free_prob(&tsp);

    return rval;
}
