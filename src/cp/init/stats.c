#include "cp/init/init.h"
#include "op-solver.h"

cp_init_stats *
cp_create_init_stats(void)
{
    cp_init_stats *stats = malloc(sizeof(cp_init_stats));

    stats->total = stats_create("CP Initialization");
    stats->write = 1;
    stats->file  = "stats.json";
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

int
cp_write_init_stats(cp_prob *cp, cp_init_env *env)
{
    int rval;
    cp_init_stats *stats = env->stats;
    cp_init_param *param = env->param;
    FILE *file           = NULL;
    struct timespec timestamp;

    printf("\n");
    printf("Writing cp_init_stats to '%s'...\n", stats->file);
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
#if HAVE_LP_SOLVER
    fprintf(file, "\"lb\": %.f, ", cp->ip->lowerboundG);
    fprintf(file, "\"ub\": %.f ", cp->ip->upperboundG);
#else
    fprintf(file, "\"lb\": %.f, ", cp->sol->val);
    fprintf(file, "\"ub\": %.f ", SOLVER_MAXDOUBLE);
#endif
    fprintf(file, "}, ");

    // Parameters
    fprintf(file, "\"param\": { ");
    fprintf(file, "\"time_limit\": %ld, ", param->time_limit);
    fprintf(file, "\"init\": %d, ", param->init);
    fprintf(file, "\"select\": %d, ", param->select);
    fprintf(file, "\"pinit\": %.2f ", param->pinit);
    fprintf(file, "}, ");

    // Stats
    fprintf(file, "\"stats\": { ");
    fprintf(file, "\"time\": %ld ", stats_get_total_time(stats->total));
    fprintf(file, "}, ");

    clock_gettime(CLOCK_REALTIME, &timestamp);
    fprintf(file, "\"timestamp\": %ld, ",
            timestamp.tv_sec * 1000 + timestamp.tv_nsec / 1000000);

    fprintf(file, "\"event\": \"%s\", ", "stats_summary");
    fprintf(file, "\"env\": \"%s\", ", "cp_init");
    fprintf(file, "\"seed\": %lu, ", __seed__);
    fprintf(file, "\"pid\": %lu ", (long unsigned)getpid());
    fprintf(file, "} \n");
    rval = 0;
done:
    if (file)
        fclose(file);
    return rval;
}
