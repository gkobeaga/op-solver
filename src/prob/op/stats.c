#include "op-solver.h"
#include "op/op.h"

op_stats *
op_create_stats(void)
{
    op_stats *stats = malloc(sizeof(op_stats));

    stats->total = stats_create("Orienteering Problem");
    return stats;
}

void
op_free_stats(op_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
