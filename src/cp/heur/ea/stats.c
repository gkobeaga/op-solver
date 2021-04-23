#include "cp/heur/ea/ea.h"
#include "op-solver.h"

cp_heur_ea_stats *
cp_create_heur_ea_stats(void)
{
    cp_heur_ea_stats *stats = malloc(sizeof(cp_heur_ea_stats));

    stats->total = stats_create("CP Heurisitc: Evolutionary Algorithm");
    stats->it    = stats_create("Iterations");
    stats->infeas_recover = stats_create("Infeas recover");
    stats->write_stats    = 0;
    stats->file           = "stats.json";
    return stats;
}

void
cp_free_heur_ea_stats(cp_heur_ea_stats **stats)
{
    if (*stats)
    {
        free((*stats)->total);
        free((*stats)->it);
        free((*stats)->infeas_recover);
        free(*stats);
        *stats = NULL;
    }
}

int
cp_write_heur_ea_stats(cp_prob *cp, cp_heur_ea_env *env)
{
    int rval;
    cp_heur_ea_stats *stats = env->stats;
    cp_heur_ea_param *param = env->param;
    FILE *file              = NULL;
    struct timespec timestamp;

    printf("\n");
    printf("Writing cp_heur_ea stats to '%s'...\n\n", stats->file);
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
    if (cp->ip->sol)
    {
        fprintf(file, "\"lb\": %.f, ", cp->ip->lowerboundG);
        fprintf(file, "\"ub\": %.f, ", cp->ip->upperboundG);
    }
    else
    {
        fprintf(file, "\"lb\": %.f, ", cp->sol->val);
        fprintf(file, "\"ub\": %.f, ", SOLVER_MAXDOUBLE);
    }
    fprintf(file, "\"cycle\": [ ");
    for (int i = 0; i < cp->sol->ns - 1; i++)
        fprintf(file, "%d, ", cp->sol->cycle[i] + 1);
    fprintf(file, "%d]", cp->sol->cycle[cp->sol->ns - 1] + 1);
    fprintf(file, "}, ");

    // Parameters
    fprintf(file, "\"param\": { ");
    fprintf(file, "\"time_limit\": %ld, ", param->time_limit);
    fprintf(file, "\"it_lim\": %d, ", param->it_lim);
    fprintf(file, "\"pop_size\": %d, ", param->pop_size);
    fprintf(file, "\"pop_stop\": %d, ", param->pop_stop);
    fprintf(file, "\"d2d\": %d, ", param->d2d);
    fprintf(file, "\"nparsel\": %d, ", param->nparsel);
    fprintf(file, "\"pmut\": %f, ", param->pmut);
    fprintf(file, "\"len_improve1\": %d, ", param->len_improve1);
    fprintf(file, "\"len_improve2\": %d ", param->len_improve2);
    fprintf(file, "}, ");

    // Stats
    fprintf(file, "\"stats\": { ");
    fprintf(file, "\"time\": %ld, ", stats_get_total_time(stats->total));
    fprintf(file, "\"it\": %d, ", stats->it->count_active);
    fprintf(file, "\"time_infeas_recover\": %ld ",
            stats_get_total_time(stats->infeas_recover));
    fprintf(file, "}, ");

    clock_gettime(CLOCK_REALTIME, &timestamp);
    fprintf(file, "\"timestamp\": %ld, ",
            timestamp.tv_sec * 1000 + timestamp.tv_nsec / 1000000);

    fprintf(file, "\"event\": \"%s\", ", "stats_summary");
    fprintf(file, "\"env\": \"%s\", ", "cp_heur_ea");
    fprintf(file, "\"seed\": %lu, ", __seed__);
    fprintf(file, "\"pid\": %lu ", (long unsigned)getpid());
    fprintf(file, "} \n");
    rval = 0;
done:
    if (file)
        fclose(file);
    return rval;
}
