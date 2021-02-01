#include "op-solver.h"
#include "op/init/init.h"

op_init_stats *
op_create_init_stats(void)
{
    op_init_stats *stats = malloc(sizeof(op_init_stats));

    stats->total = stats_create("OP Initialization");
    return stats;
}

void
op_free_init_stats(op_init_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
