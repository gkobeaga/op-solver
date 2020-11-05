#include "op-solver.h"
#include "tsp/tsp.h"

tsp_stats *
tsp_create_stats(void)
{
    tsp_stats *stats = malloc(sizeof(tsp_stats));

    stats->total = stats_create("Integer Problem");
    return stats;
}

void
tsp_free_stats(tsp_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
