#include "cp/heur/ea/ea.h"
#include "op-solver.h"

cp_heur_ea_stats *
cp_create_heur_ea_stats(void)
{
    cp_heur_ea_stats *stats = malloc(sizeof(cp_heur_ea_stats));

    stats->total = stats_create("CP Heurisitc: Evolutionary Algorithm");
    stats->it    = stats_create("Iterations");
    stats->infeas_recover = stats_create("Infeas recover");
    return stats;
}

void
cp_free_heur_ea_stats(cp_heur_ea_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free((*stats)->it);
        free((*stats)->infeas_recover);
        free(*stats);
        *stats = NULL;
    }
}
