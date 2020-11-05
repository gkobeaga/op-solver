#include "op-solver.h"

int
rng_bernoulli(double p)
{
    assert(0.0 <= p && p <= 1.0);

    if (drand48() < p)
    {
        return 1;
    }
    {
        return 0;
    }
}

static void
copy(void *dest, size_t i, void *src, size_t j, size_t size)
{
    char *a  = size * i + (char *)dest;
    char *b  = size * j + (char *)src;
    size_t s = size;

    do
    {
        *a++ = *b++;
    } while (--s > 0);
}

void
rng_choose(void *candidates, int ncand, void *selected, int nsel, size_t size)
{
    int i, j = 0;

    assert(nsel <= ncand);

    for (i = 0; i < ncand && j < nsel; i++)
    {
        if ((rand() % (ncand - i)) < nsel - j)
        {
            copy(selected, j, candidates, i, size);
            j++;
        }
    }
}

static int
compare(const void *a, const void *b)
{
    if (*(double *)a > *(double *)b)
        return 1;
    else if (*(double *)a < *(double *)b)
        return -1;
    return 0;
}

void
rng_choose_p(void *candidates, double *p, int ncand, void *selected, int nsel,
             size_t size)
{
    double cum;
    int i, j;
    double *x = NULL;

    for (i = 0, cum = 0; i < ncand; i++)
    {
        cum += p[i];
    }

    assert(cum >= 1.0 - SOLVER_ZEROPLUS);
    assert(cum <= 1.0 + SOLVER_ZEROPLUS);

    x = malloc(nsel * sizeof(double));

    for (j = 0; j < nsel; j++)
    {
        x[j] = drand48();
    }

    qsort(x, nsel, sizeof(double), compare);

    cum = 0.0;
    i = j = 0;
    do
    {
        if (x[j] <= cum + p[i])
        {
            copy(selected, j, candidates, i, size);
            j++;
        }
        else
        {
            cum += p[i];
            i++;
        }

    } while (i < ncand && j < nsel);

    if (x)
        free(x);
}
