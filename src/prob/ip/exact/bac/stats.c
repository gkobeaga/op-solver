#include "ip/exact/bac/bac.h"
#include "op-solver.h"

ip_exact_bac_stats *
ip_create_exact_bac_stats(void)
{
    ip_exact_bac_stats *stats = malloc(sizeof(ip_exact_bac_stats));

    stats->total                = stats_create("IP Exact Approach");
    stats->branch_node          = stats_create("Branch node");
    stats->pricing_loop         = stats_create("Pricing Loop");
    stats->find_branch          = stats_create("Find branch");
    stats->exec_branch          = stats_create("Exec branch");
    stats->exec_unbranch        = stats_create("Exec unbranch");
    stats->verify_branch_infeas = stats_create("Verify Infeasible Branch");
    stats->verify_branch_prune  = stats_create("Verify Branch Prune");
    stats->recover_infeas       = stats_create("Recover Infeasible");

    stats->sep_loop          = stats_create("Separation Loop");
    stats->sparse_edge_check = stats_create("Sparse edge check");
    stats->full_edge_check   = stats_create("Full edge check");
    stats->addcuts           = stats_create("Add cuts");
    stats->agecuts           = stats_create("Age cuts");
    stats->ageedges          = stats_create("Age edges");
    stats->addbad            = stats_create("Add edges");
    stats->xheur             = stats_create("X-heuristic");
    stats->misc              = stats_create("Miscellaneous");
    return stats;
}

void
ip_free_exact_bac_stats(ip_exact_bac_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);

        free((*stats)->pricing_loop);
        free((*stats)->branch_node);
        free((*stats)->find_branch);
        free((*stats)->exec_branch);
        free((*stats)->exec_unbranch);
        free((*stats)->verify_branch_infeas);
        free((*stats)->verify_branch_prune);
        free((*stats)->recover_infeas);

        free((*stats)->sep_loop);
        free((*stats)->sparse_edge_check);
        free((*stats)->full_edge_check);
        free((*stats)->addcuts);
        free((*stats)->agecuts);
        free((*stats)->ageedges);
        free((*stats)->addbad);
        free((*stats)->xheur);
        free((*stats)->misc);
        free(*stats);
        *stats = NULL;
    }
}

// TODO
int
ip_write_exact_bac_stats(ip_exact_bac_env *env)
{
    return 0;
}
