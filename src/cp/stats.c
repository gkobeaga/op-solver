#include "cp/cp.h"
#include "op-solver.h"

cp_stats *
cp_create_stats(void)
{
    cp_stats *stats = malloc(sizeof(cp_stats));

    stats->total = stats_create("Cycle Problem");
    return stats;
}

void
cp_free_stats(cp_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
