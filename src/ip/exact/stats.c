#include "ip/exact/exact.h"
#include "op-solver.h"

ip_exact_stats *
ip_create_exact_stats(void)
{
    ip_exact_stats *stats = malloc(sizeof(ip_exact_stats));

    stats->total = stats_create("IP Exact Approach");
    return stats;
}

void
ip_free_exact_stats(ip_exact_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
