#include "op-solver.h"

struct stats_item *
stats_create(const char *name)
{
    struct stats_item *item = malloc(sizeof(struct stats_item));
    item->active            = 0;
    item->total_time        = 0.0;
    item->count_active      = 0;
    item->count_success     = 0;
    item->count_total       = 0;
    item->last_success      = 0;
    item->last_total        = 0;
    if (name == (char *)NULL || name[0] == '\0')
        strncpy(item->name, "ANONYMOUS", sizeof(item->name) - 1);
    else
        strncpy(item->name, name, sizeof(item->name) - 1);
    item->name[sizeof(item->name) - 1] = '\0';
    return item;
}

int
stats_start(struct stats_item *item)
{
    if (item->active != 0)
    {
        fprintf(stderr, "Warning: restarting running item %s\n", item->name);
        return 1;
    }
    clock_gettime(CLOCK_REALTIME, &item->current_start);
    item->active = 1;

    return 0;
}

int
stats_suspend(struct stats_item *item)
{
    if (item->active == 0)
    {
        fprintf(stderr, "Warning: suspended non-running item %s\n", item->name);
        return 1;
    }

    clock_gettime(CLOCK_REALTIME, &item->current_end);

    clock_gettime(CLOCK_REALTIME, &item->last_end);
    item->total_time +=
    (long)((item->current_end.tv_sec - item->current_start.tv_sec) * 1000000 +
           (item->current_end.tv_nsec - item->current_start.tv_nsec) / 1000);
    item->active = 0;

    return 0;
}

int
stats_resume(struct stats_item *item)
{
    if (item->active != 0)
    {
        fprintf(stderr, "Warning: resuming running item %s\n", item->name);
        return 1;
    }
    clock_gettime(CLOCK_REALTIME, &item->last_end);
    item->active = 1;

    return 0;
}

int
stats_stop(struct stats_item *item, int count)
{
    if (item->active == 0)
    {
        fprintf(stderr, "Warning: stopping non-running item %s\n", item->name);
        return 1;
    }

    clock_gettime(CLOCK_REALTIME, &item->last_end);

    item->total_time +=
    (long)((item->last_end.tv_sec - item->current_start.tv_sec) * 1000000 +
           (item->last_end.tv_nsec - item->current_start.tv_nsec) / 1000);
    item->active = 0;
    item->count_active++;
    item->count_total += count;
    if (count)
    {
        item->count_success++;
        item->last_success = 1;
        item->last_total   = count;
    }
    else
    {
        item->last_success = 0;
        item->last_total   = 0;
    }

    return 0;
}

int
stats_stop_if_active(struct stats_item *item, int count)
{
    if (item->active == 0)
    {
        return 0;
    }

    clock_gettime(CLOCK_REALTIME, &item->last_end);

    item->total_time +=
    (long)((item->last_end.tv_sec - item->current_start.tv_sec) * 1000000 +
           (item->last_end.tv_nsec - item->current_start.tv_nsec) / 1000);
    item->active = 0;
    item->count_active++;
    item->count_total += count;
    if (count)
    {
        item->count_success++;
        item->last_success = 1;
        item->last_total   = count;
    }
    else
    {
        item->last_success = 0;
        item->last_total   = 0;
    }

    return 0;
}

long
stats_get_current_time(struct stats_item *item)
{
    struct timespec current;
    clock_gettime(CLOCK_REALTIME, &current);
    long z = (long)((current.tv_sec - item->current_start.tv_sec) * 1000000 +
                    (current.tv_nsec - item->current_start.tv_nsec) / 1000);
    return z / 1000;
}

long
stats_get_total_time(struct stats_item *item)
{
    struct timespec current;

    if (item->active != 0)
    {
        clock_gettime(CLOCK_REALTIME, &current);
        long z =
        (long)((current.tv_sec - item->current_start.tv_sec) * 1000000 +
               (current.tv_nsec - item->current_start.tv_nsec) / 1000);
        return z / 1000;
    }
    else
    {
        return item->total_time / 1000;
    }
}

int
stats_print(struct stats_item *item)
{
    struct timespec current;

    printf("Stats for %-34.34s: \n", item->name);
    if (item->active != 0)
    {
        clock_gettime(CLOCK_REALTIME, &current);
        long z =
        (long)((current.tv_sec - item->current_start.tv_sec) * 1000000 +
               (current.tv_nsec - item->current_start.tv_nsec) / 1000);
        printf("  Time: active %.2ld ms, total %.2ld ms\n", z,
               z + item->total_time);
    }
    else
    {
        printf("  Count: success %d/%d total %d/%d active %d calls\n",
               item->last_success, item->count_success, item->last_total,
               item->count_total, item->count_active);
        printf("  Time:  last %.2ld ms, cumulative %.2ld ms\n", item->last_time,
               item->total_time);
    }

    return 0;
}
