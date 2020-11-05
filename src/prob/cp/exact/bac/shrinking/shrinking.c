#include "../bac.h"
#include "op-solver.h"

int
cp_shrink_exact_bac_graph(cp_prob *cp, cp_exact_bac_env *bac_env,
                          solver_graph *graph, graph_vertex *qstart,
                          graph_clique_repo *repo)
{
    int rval = 0;

    switch (bac_env->param->srk_rule)
    {
    case CP_SRK_C1:
        cp_shrink_exact_bac_graph_c1(cp, bac_env, graph, qstart, repo);
        break;
    case CP_SRK_C1C2:
        cp_shrink_exact_bac_graph_c1c2(cp, bac_env, graph, qstart, repo);
        break;
    case CP_SRK_C1C2C3:
        cp_shrink_exact_bac_graph_c1c2c3(cp, bac_env, graph, qstart, repo);
        break;
    case CP_SRK_S1:
        cp_shrink_exact_bac_graph_s1(cp, bac_env, graph, qstart, repo);
        break;
    default:
        rval = 1;
        fprintf(stderr, "Invalid shrinking rule in cp_shrink_exact_bac\n");
        break;
    }

    return rval;
}
