#include "cp/exact/exact.h"
#include "op-solver.h"

cp_exact_stats *
cp_create_exact_stats(void)
{
    cp_exact_stats *stats = malloc(sizeof(cp_exact_stats));

    stats->total = stats_create("CP Exact Approach");
    return stats;
}

void
cp_free_exact_stats(cp_exact_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
