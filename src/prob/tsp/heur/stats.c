#include "op-solver.h"
#include "tsp/heur/heur.h"

tsp_heur_stats *
tsp_create_heur_stats(void)
{
    tsp_heur_stats *stats = malloc(sizeof(tsp_heur_stats));

    stats->total = stats_create("TSP Heuristic Approach");
    return stats;
}

void
tsp_free_heur_stats(tsp_heur_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
