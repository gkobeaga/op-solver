#include "cp/heur/heur.h"
#include "op-solver.h"

cp_heur_stats *
cp_create_heur_stats(void)
{
    cp_heur_stats *stats = malloc(sizeof(cp_heur_stats));

    stats->total = stats_create("CP Heuristic Approach");
    return stats;
}

void
cp_free_heur_stats(cp_heur_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
