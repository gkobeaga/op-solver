#include "op-solver.h"
#include "cp/exact/bac/bac.h"

cp_exact_bac_stats *
cp_create_exact_bac_stats(void)
{
    cp_exact_bac_stats *stats = malloc(sizeof(cp_exact_bac_stats));

    stats->total             = stats_create("CP Branch-and-Cut");
    stats->sep_loop          = stats_create("Separation Loop");
    stats->sep_loop_it       = stats_create("Separation loop iteration");
    stats->sep_loop_inner    = stats_create("Inner separation loop");
    stats->sep_loop_inner_it = stats_create("Inner separation loop iteration");
    stats->sep_loop_middle   = stats_create("Middle separation loop");
    stats->sep_loop_middle_it =
    stats_create("Middle separation loop iteration");
    stats->sep_loop_outer     = stats_create("Outer separation loop");
    stats->sep_loop_outer_it  = stats_create("Outer separation loop iteration");
    stats->sep_logical        = stats_create("Logical");
    stats->sep_sec_comps      = stats_create("Component SEC");
    stats->sep_sec_exact      = stats_create("Exact SEC");
    stats->sep_connect_mincut = stats_create("Connectiviy");
    stats->sep_blossom_mst    = stats_create("MST blossom");
    stats->sep_blossom_fast   = stats_create("Fast blossom");
    stats->sep_blossom_ghfast = stats_create("Fast GH-blossom");
    stats->sep_cover_edge     = stats_create("Edge cover");
    stats->sep_cover_vertex   = stats_create("Vertex cover");
    stats->sep_cover_cycle    = stats_create("Cycle cover");
    stats->sep_path           = stats_create("Path");
    stats->age_cuts           = stats_create("Age cuts");
    stats->age_vars           = stats_create("Age edges");
    stats->add_cuts           = stats_create("Add cuts");
    stats->add_vars           = stats_create("Add variables");
    stats->xheur_branch       = stats_create("Separation x-heuristic");
    stats->xheur_sep          = stats_create("Branch x-heuristic");
    stats->lp_opt             = stats_create("LP opt");
    stats->misc               = stats_create("Miscellaneous");
    stats->file               = "bac-stats.json";
    return stats;
}

void
cp_free_exact_bac_stats(cp_exact_bac_stats **stats)
{
    if (*stats)
    {

        free((*stats)->total);
        free((*stats)->sep_loop);
        free((*stats)->sep_loop_it);
        free((*stats)->sep_loop_inner);
        free((*stats)->sep_loop_inner_it);
        free((*stats)->sep_loop_middle);
        free((*stats)->sep_loop_middle_it);
        free((*stats)->sep_loop_outer);
        free((*stats)->sep_loop_outer_it);
        free((*stats)->sep_logical);
        free((*stats)->sep_sec_comps);
        free((*stats)->sep_sec_exact);
        free((*stats)->sep_blossom_fast);
        free((*stats)->sep_blossom_ghfast);
        free((*stats)->sep_blossom_mst);
        free((*stats)->sep_connect_mincut);
        free((*stats)->sep_cover_edge);
        free((*stats)->sep_cover_vertex);
        free((*stats)->sep_cover_cycle);
        free((*stats)->sep_path);
        free((*stats)->age_cuts);
        free((*stats)->age_vars);
        free((*stats)->add_cuts);
        free((*stats)->add_vars);
        free((*stats)->xheur_branch);
        free((*stats)->xheur_sep);
        free((*stats)->lp_opt);
        free((*stats)->misc);
        free(*stats);
        *stats = NULL;
    }
}

int
cp_write_exact_bac_stats(cp_prob *cp, cp_exact_bac_env *env)
{
    int rval;
    cp_exact_bac_stats *stats = env->stats;
    cp_exact_bac_param *param = env->param;
    FILE *file                = NULL;
    struct timespec timestamp;

    printf("\n");
    printf("Writing excution stats to '%s'...\n", stats->file);
    if (cp == NULL)
    {
        printf("No CP problem to write.\n");
    }

    file = fopen(stats->file, "a");
    check_null(file, "Unable to create file", done);

    fprintf(file, "{ ");

    // Problem
    fprintf(file, "\"prob\": { ");
    fprintf(file, "\"name\": \"%s\", ", cp->data->name);
    fprintf(file, "\"n\": %d, ", cp->n);
    fprintf(file, "\"d0\": %.0f ", cp->data->cap);
    fprintf(file, "}, ");

    // Solution
    fprintf(file, "\"sol\": { ");
    fprintf(file, "\"val\": %.0f, ", cp->sol->val);
    fprintf(file, "\"cap\": %.0f, ", cp->sol->cap);
    fprintf(file, "\"sol_ns\": %d, ", cp->sol->ns);
    fprintf(file, "\"lb\": %.f, ", cp->ip->lowerboundG);
    fprintf(file, "\"ub\": %.f ", cp->ip->upperboundG);
    fprintf(file, "}, ");

    // Parameters
    fprintf(file, "\"param\": { ");
    fprintf(file, "\"sep_logical\": %d, ", param->sep_logical);
    fprintf(file, "\"sep_sec_comps\": %d, ", param->sep_sec_comps);
    fprintf(file, "\"sep_sec_exact\": %d, ", param->sep_sec_exact);
    fprintf(file, "\"sep_sec_cc_2\": %d, ", param->sep_sec_cc_2);
    fprintf(file, "\"sep_sec_cc_extra\": %d, ", param->sep_sec_cc_extra);
    fprintf(file, "\"sep_blossom_fst\": %d, ", param->sep_blossom_mst);
    fprintf(file, "\"sep_blossom_eph\": %d, ", param->sep_blossom_fast);
    fprintf(file, "\"sep_blossom_egh\": %d, ", param->sep_blossom_ghfast);
    fprintf(file, "\"sep_cover_edge\": %d, ", param->sep_cover_edge);
    fprintf(file, "\"sep_cover_vertex\": %d, ", param->sep_cover_vertex);
    fprintf(file, "\"sep_cover_cycle\": %d, ", param->sep_cover_cycle);
    fprintf(file, "\"sep_path\": %d, ", param->sep_path);
    fprintf(file, "\"sep_loop\": %d, ", param->sep_loop);
    fprintf(file, "\"sep_srk_rule\": %d, ", param->srk_rule);
    fprintf(file, "\"sep_srk_s2\": %d, ", param->srk_s2);
    fprintf(file, "\"sep_srk_s3\": %d, ", param->srk_s3);
    fprintf(file, "\"sep_srk_extra\": %d, ", param->srk_extra);
    fprintf(file, "\"xheur_vph\": %d, ", param->xheur_vph);
    fprintf(file, "\"xheur_vph_meta\": %d ", param->xheur_vph_meta);
    fprintf(file, "}, ");

    // Stats
    fprintf(file, "\"stats\": { ");
    fprintf(file, "\"time\": %ld, ", stats_get_total_time(stats->total));
    fprintf(file, "\"active_sep_logical\": %d, ",
            stats->sep_logical->count_active);
    fprintf(file, "\"success_sep_logical\": %d, ",
            stats->sep_logical->count_success);
    fprintf(file, "\"total_sep_logical\": %d, ",
            stats->sep_logical->count_total);

    fprintf(file, "\"active_sep_sec_comps\": %d, ",
            stats->sep_sec_comps->count_active);
    fprintf(file, "\"success_sep_sec_comps\": %d, ",
            stats->sep_sec_comps->count_success);
    fprintf(file, "\"total_sep_sec_comps\": %d, ",
            stats->sep_sec_comps->count_total);
    fprintf(file, "\"time_sep_sec_comps\": %ld, ",
            stats_get_total_time(stats->sep_sec_comps));

    fprintf(file, "\"active_sep_sec_exact\": %d, ",
            stats->sep_sec_exact->count_active);
    fprintf(file, "\"success_sep_sec_exact\": %d, ",
            stats->sep_sec_exact->count_success);
    fprintf(file, "\"total_sep_sec_exact\": %d, ",
            stats->sep_sec_exact->count_total);
    fprintf(file, "\"time_sep_sec_exact\": %ld, ",
            stats_get_total_time(stats->sep_sec_exact));

    fprintf(file, "\"active_sep_blossom_fast\": %d, ",
            stats->sep_blossom_fast->count_active);
    fprintf(file, "\"success_sep_blossom_fast\": %d, ",
            stats->sep_blossom_fast->count_success);
    fprintf(file, "\"total_sep_blossom_fast\": %d, ",
            stats->sep_blossom_fast->count_total);
    fprintf(file, "\"time_sep_blossom_fast\": %ld, ",
            stats_get_total_time(stats->sep_blossom_fast));

    fprintf(file, "\"active_sep_blossom_ghfast\": %d, ",
            stats->sep_blossom_ghfast->count_active);
    fprintf(file, "\"success_sep_blossom_ghfast\": %d, ",
            stats->sep_blossom_ghfast->count_success);
    fprintf(file, "\"total_sep_blossom_ghfast\": %d, ",
            stats->sep_blossom_ghfast->count_total);
    fprintf(file, "\"time_sep_blossom_fast\": %ld, ",
            stats_get_total_time(stats->sep_blossom_ghfast));

    fprintf(file, "\"active_sep_blossom_mst\": %d, ",
            stats->sep_blossom_mst->count_active);
    fprintf(file, "\"success_sep_blossom_mst\": %d, ",
            stats->sep_blossom_mst->count_success);
    fprintf(file, "\"total_sep_blossom_mst\": %d, ",
            stats->sep_blossom_mst->count_total);
    fprintf(file, "\"time_sep_blossom_mst\": %ld, ",
            stats_get_total_time(stats->sep_blossom_mst));

    fprintf(file, "\"active_sep_cover_edge\": %d, ",
            stats->sep_cover_edge->count_active);
    fprintf(file, "\"success_sep_cover_edge\": %d, ",
            stats->sep_cover_edge->count_success);
    fprintf(file, "\"total_sep_cover_edge\": %d, ",
            stats->sep_cover_edge->count_total);
    fprintf(file, "\"time_sep_cover_edge\": %ld, ",
            stats_get_total_time(stats->sep_cover_edge));

    fprintf(file, "\"active_sep_cover_cycle\": %d, ",
            stats->sep_cover_cycle->count_active);
    fprintf(file, "\"success_sep_cover_cycle\": %d, ",
            stats->sep_cover_cycle->count_success);
    fprintf(file, "\"total_sep_cover_cycle\": %d, ",
            stats->sep_cover_cycle->count_total);
    fprintf(file, "\"time_sep_cover_cycle\": %ld, ",
            stats_get_total_time(stats->sep_cover_cycle));

    fprintf(file, "\"active_sep_cover_vertex\": %d, ",
            stats->sep_cover_vertex->count_active);
    fprintf(file, "\"success_sep_cover_vertex\": %d, ",
            stats->sep_cover_vertex->count_success);
    fprintf(file, "\"total_sep_cover_vertex\": %d, ",
            stats->sep_cover_edge->count_total);
    fprintf(file, "\"time_sep_cover_vertex\": %ld, ",
            stats_get_total_time(stats->sep_cover_vertex));

    fprintf(file, "\"active_sep_path\": %d, ", stats->sep_path->count_active);
    fprintf(file, "\"success_sep_path\": %d, ", stats->sep_path->count_success);
    fprintf(file, "\"total_sep_path\": %d, ", stats->sep_path->count_total);
    fprintf(file, "\"time_sep_path\": %ld, ",
            stats_get_total_time(stats->sep_path));

    fprintf(file, "\"time_sep_sep_loop\": %ld, ",
            stats_get_total_time(stats->sep_loop));
    fprintf(file, "\"time_sep_sep_loop_it\": %ld, ",
            stats_get_total_time(stats->sep_loop_it));
    fprintf(file, "\"time_sep_sep_loop_inner\": %ld, ",
            stats_get_total_time(stats->sep_loop_inner));
    fprintf(file, "\"time_sep_sep_loop_inner_it\": %ld, ",
            stats_get_total_time(stats->sep_loop_inner_it));
    fprintf(file, "\"time_sep_sep_loop_middle\": %ld, ",
            stats_get_total_time(stats->sep_loop_middle));
    fprintf(file, "\"time_sep_sep_loop_middle_it\": %ld, ",
            stats_get_total_time(stats->sep_loop_middle_it));
    fprintf(file, "\"time_sep_sep_loop_outer\": %ld, ",
            stats_get_total_time(stats->sep_loop_outer));
    fprintf(file, "\"time_sep_sep_loop_outer_it\": %ld, ",
            stats_get_total_time(stats->sep_loop_outer_it));
    fprintf(file, "\"time_age_cut\": %ld, ",
            stats_get_total_time(stats->age_cuts));
    fprintf(file, "\"time_age_vars\": %ld, ",
            stats_get_total_time(stats->age_vars));
    fprintf(file, "\"time_add_vars\": %ld, ",
            stats_get_total_time(stats->add_vars));
    fprintf(file, "\"time_add_cuts\": %ld, ",
            stats_get_total_time(stats->add_cuts));
    fprintf(file, "\"time_xheur_branch\": %ld, ",
            stats_get_total_time(stats->xheur_branch));
    fprintf(file, "\"time_xheur_sep\": %ld ",
            stats_get_total_time(stats->xheur_sep));
    fprintf(file, "}, ");

    clock_gettime(CLOCK_REALTIME, &timestamp);
    fprintf(file, "\"timestamp\": %ld, ",
            timestamp.tv_sec * 1000 + timestamp.tv_nsec / 1000000);

    fprintf(file, "\"event\": \"%s\", ", "stats_summary");
    fprintf(file, "\"env\": \"%s\", ", "cp_exact_bac");
    fprintf(file, "\"seed\": %lu, ", __seed__);
    fprintf(file, "\"pid\": %lu ", (long unsigned)getpid());
    fprintf(file, "} \n");
    rval = 0;
done:
    if (file)
        fclose(file);
    return rval;
}
