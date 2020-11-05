#include "op-solver.h"
#include "tsp/heur/linkern/linkern.h"

tsp_heur_linkern_stats *
tsp_create_heur_linkern_stats(void)
{
    tsp_heur_linkern_stats *stats = malloc(sizeof(tsp_heur_linkern_stats));

    stats->total = stats_create("TSP Lin-Kernighan Approach");
    stats->steps = stats_create("Lin-Kernighan steps");
    return stats;
}

void
tsp_free_heur_linkern_stats(tsp_heur_linkern_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free((*stats)->steps);
        free(*stats);
        *stats = NULL;
    }
}
