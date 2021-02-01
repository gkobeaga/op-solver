#include "ip/ip.h"
#include "op-solver.h"

ip_stats *
ip_create_stats(void)
{
    ip_stats *stats = malloc(sizeof(ip_stats));

    stats->total = stats_create("Integer Problem");
    return stats;
}

void
ip_free_stats(ip_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free(*stats);
        *stats = NULL;
    }
}
