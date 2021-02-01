#include "cp/cp.h"
#include "op-solver.h"

static void
cp_create_pop_work(cp_pop *pop, int size)
{
    int i;
    pop->size     = size;
    pop->sol      = malloc(size * sizeof(cp_sol *));
    pop->rankperm = malloc(size * sizeof(int));
    for (i = 0; i < size; i++)
    {
        pop->sol[i]      = NULL;
        pop->rankperm[i] = i;
    }
    pop->mean_val  = 0.0;
    pop->best_val  = 0.0;
    pop->best_ind  = 0;
    pop->q25_ind   = 0;
    pop->q25_val   = 0.0;
    pop->q50_ind   = 0;
    pop->q50_val   = 0.0;
    pop->q75_ind   = 0;
    pop->q75_val   = 0.0;
    pop->stop_val  = 0.0;
    pop->stop_ind  = 0;
    pop->worst_val = 0.0;
    pop->worst_ind = 0;
    pop->parent    = NULL;
}

cp_pop *
cp_create_pop(cp_prob *cp, int size)
{
    int i;
    cp_pop *pop = malloc(sizeof(cp_pop));
    cp_create_pop_work(pop, size);
    for (i = 0; i < size; i++)
    {
        pop->sol[i] = cp_create_sol(cp);
    }
    return pop;
}

void
cp_set_pop_sol(cp_pop *pop, cp_sol *sol, int pos)
{
    cp_sol *worst = pop->sol[pop->worst_ind];
    cp_copy_sol(sol, worst);
    return;
}

static int
sort_sols(const void *ss, const void *tt)
{
    cp_sol *s = *(cp_sol **)ss, *t = *(cp_sol **)tt;

    if (s->val < t->val)
        return 1;
    else
        return -1;
    return 0;
}

int
cp_update_pop(cp_pop *pop)
{
    int rval = 0;
    int i, qstep;
    pop->mean_val = 0.0;

    for (i = 0; i < pop->size; i++)
    {
        cp_sol *sol      = pop->sol[i];
        pop->rankperm[i] = i;
        pop->mean_val += sol->val / pop->size;
    }

    qsort(pop->sol, pop->size, sizeof(cp_sol *), sort_sols);

    pop->best_ind = 0;
    pop->best_val = pop->sol[0]->val;
    qstep         = floor(pop->size / 4.0);
    pop->q25_ind  = qstep;
    pop->q25_val  = pop->sol[pop->q25_ind]->val;
    pop->q50_ind  = 2 * qstep;
    pop->q50_val  = pop->sol[pop->q50_ind]->val;
    pop->q75_ind  = 3 * qstep;
    pop->q75_val  = pop->sol[pop->q75_ind]->val;

    pop->worst_ind = pop->size - 1;
    pop->worst_val = pop->sol[pop->worst_ind]->val;

    return rval;
}

void
cp_delete_pop(cp_pop *pop)
{
    int i;

    if (pop->parent)
        free(pop->parent);
    for (i = 0; i < pop->size; i++)
    {
        cp_free_sol(&pop->sol[i]);
    }
    free(pop->sol);
    free(pop->rankperm);
}

void
cp_erase_pop(cp_pop *pop)
{
    int size;
    size = pop->size;
    cp_delete_pop(pop);
    cp_create_pop_work(pop, size);
    return;
}

void
cp_free_pop(cp_pop **pop)
{
    cp_delete_pop(*pop);
    free(*pop);
    *pop = NULL;
    return;
}
