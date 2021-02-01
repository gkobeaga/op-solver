#include "kp/kp.h"
#include "op-solver.h"

kp_stats *
kp_create_stats(void)
{
    kp_stats *stats = malloc(sizeof(kp_stats));

    stats->total       = stats_create("Integer Problem");
    stats->init_sol    = stats_create("Initial Solution");
    stats->exact       = stats_create("Exact Algorithm");
    stats->branch_node = stats_create("Branch Node");
    return stats;
}

void
kp_free_stats(kp_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free((*stats)->init_sol);
        free((*stats)->exact);
        free((*stats)->branch_node);
        free(*stats);
        *stats = NULL;
    }
}
