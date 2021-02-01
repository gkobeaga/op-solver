#include "op-solver.h"
#include "ip/ip.h"
#include "ip/exact/bac/bac.h"

void
ip_print_branch_history(ip_exact_bac_env *bac_env)
{
    int j;
    ip_branch *branch;
    printf("Branch History\n");
    if (bac_env->history_depth == 0)
    {
        printf("    Root Node\n");
    }
    else
    {
        for (j = 0; j < bac_env->history_depth; j++)
        {
            printf("    ");
            branch = bac_env->history[j];
            printf("Edge (%d,%d) set to %d\n", branch->edge->tail->i,
                   branch->edge->head->i, branch->rhs);
        }
    }
}
