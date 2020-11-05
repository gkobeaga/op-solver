#include "ip/ip.h"
#include "op-solver.h"

static void
ip_create_prob_work(ip_prob *prob)
{
    prob->graph = NULL;
    prob->cuts  = NULL;
    prob->lp    = lp_create_prob();
    prob->sense = SOLVER_OPT_SENSE_MAX;
    prob->sol   = NULL;

    prob->upperboundG = (int)SOLVER_MAXDOUBLE;
    prob->lowerboundG = (int)-SOLVER_MAXDOUBLE;
    prob->infeasible  = 0;
}

ip_prob *
ip_create_prob(void)
{
    ip_prob *prob = malloc(sizeof(ip_prob));
    ip_create_prob_work(prob);
    return prob;
}

void
ip_free_prob(ip_prob **prob)
{
    if (!(*prob))
        return;

    graph_free(&(*prob)->graph);

    lp_free_prob(&(*prob)->lp);

    free((*prob)->cuts);

    ip_free_sol(&(*prob)->sol);

    graph_free(&(*prob)->graph);

    free(*prob);
}
