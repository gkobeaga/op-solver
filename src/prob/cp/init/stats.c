#include "cp/init/init.h"
#include "op-solver.h"

cp_init_stats *
cp_create_init_stats(void)
{
    cp_init_stats *stats = malloc(sizeof(cp_init_stats));

    stats->total = stats_create("CP Initialization");
    return stats;
}

void
cp_free_init_stats(cp_init_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
