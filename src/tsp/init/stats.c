#include "op-solver.h"
#include "tsp/init/init.h"

tsp_init_stats *
tsp_create_init_stats(void)
{
    tsp_init_stats *stats = malloc(sizeof(tsp_init_stats));

    stats->total = stats_create("TSP Initialization");
    return stats;
}

void
tsp_free_init_stats(tsp_init_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
