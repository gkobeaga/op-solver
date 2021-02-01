#include "lp/lp.h"

lp_stats *
lp_create_stats(void)
{
    lp_stats *stats = malloc(sizeof(lp_stats));
    stats->total    = stats_create("Linear Problem");
    return stats;
}

void
lp_free_stats(lp_stats **stats)
{
    if (*stats)
    {
        free(*stats);
        *stats = NULL;
    }
}
