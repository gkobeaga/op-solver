#include "cp/cp.h"
#include "cp/exact/exact.h"
#include "cp/exact/bac/bac.h"

cp_exact_bac_param *
cp_create_exact_bac_param(void)
{
    cp_exact_bac_param *param   = malloc(sizeof(cp_exact_bac_param));
    param->time_limit           = 5 * 60 * 60 * 1000;
    param->sep_logical          = 1;
    param->sep_sec_comps        = 1;
    param->sep_sec_exact        = CP_SEC_HONG_BOTH;
    param->sep_sec_cc_2         = 0;
    param->sep_sec_cc_extra     = 1;
    param->sep_blossom_mst      = 0;
    param->sep_blossom_fast     = 1;
    param->sep_blossom_ghfast   = 1;
    param->sep_cover_edge       = 1;
    param->sep_cover_vertex     = 0;
    param->sep_cover_cycle      = 1;
    param->sep_path             = 1;
    param->sep_loop             = CP_SEP_LOOP_THREE_LEVEL;
    param->srk_rule             = CP_SRK_S1;
    param->srk_s2               = 0;
    param->srk_s3               = 1;
    param->srk_extra            = 1;
    param->sec_max_vout         = SOLVER_MAXINT;
    param->sec_max_vin          = SOLVER_MAXINT;
    param->sec_max_viol         = 0;
    param->sec_max_cut          = 2000;
    param->sec_max_cut_x_clique = 50;
    param->xheur_vph            = 1;
    param->xheur_vph_meta       = 1;
    param->pruning_tol          = 1.0 - SOLVER_PRICE_MAXPENALTY;
    return param;
}

void
cp_free_exact_bac_param(cp_exact_bac_param **param)
{
    if (*param)
    {
        free(*param);
        *param = NULL;
    }
}
