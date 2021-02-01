#include "op-solver.h"
#include "tsp/heur/heur.h"

int
tsp_opt_heur_twoopt(tsp_prob *tsp, tsp_heur_env *heur_env, tsp_sol *sol)
{
    solver_data *data = tsp->data;
    if (data_is_norm_type(data, SOLVER_DATA_TYPE_EUCLIDEAN))
    {
        if (heur_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("No 2-opt\n");
    }
    else if (data_is_norm_type(data, SOLVER_DATA_TYPE_HILBERT))
    {
        if (heur_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("No 2-opt\n");
    }
    else if (data_is_norm_type(data, SOLVER_DATA_TYPE_BANACH))
    {
        if (heur_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("No 2-opt\n");
    }
    return 0;
}

int
tsp_opt_heur_twoopt5(tsp_prob *tsp, tsp_heur_env *heur_env, tsp_sol *sol)
{
    solver_data *data = tsp->data;

    if (data_is_norm_type(data, SOLVER_DATA_TYPE_EUCLIDEAN))
    {
        if (heur_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("No 2-opt\n");
    }
    else if (data_is_norm_type(data, SOLVER_DATA_TYPE_HILBERT))
    {
        if (heur_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("No 2-opt\n");
    }
    else if (data_is_norm_type(data, SOLVER_DATA_TYPE_BANACH))
    {
        if (heur_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("No 2-opt\n");
    }
    return 0;
}

int
tsp_opt_heur_threeopt(tsp_prob *tsp, tsp_heur_env *heur_env, tsp_sol *sol)
{
    solver_data *data = tsp->data;
    if (data_is_norm_type(data, SOLVER_DATA_TYPE_EUCLIDEAN))
    {
        if (heur_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("No 3-opt\n");
    }
    else if (data_is_norm_type(data, SOLVER_DATA_TYPE_HILBERT))
    {
        if (heur_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("No 3-opt\n");
    }
    else if (data_is_norm_type(data, SOLVER_DATA_TYPE_BANACH))
    {
        if (heur_env->verbosity >= SOLVER_VERBOSITY_INFO)
            printf("No 3-opt\n");
    }
    else
        printf("No valid heuristic\n");
    return 0;
}
