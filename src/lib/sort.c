#include "op-solver.h"

/* NSAMPLES should be odd */
#define NSAMPLES 3
#define SORTSIZE 20

static void
select_split(int *arr, int n, double v, int *start, int *end, double *coord),
select_sort(int *arr, int n, double *coord),
select_sort_dsample(double *samp, int n);

void
sort_partial(int *arr, int lo, int hi, int mid, double *coord)
{
    double samplevals[NSAMPLES];
    int i;
    int start, end;
    int n;

    arr += lo;
    n = hi - lo + 1;
    mid -= lo;

    while (n > SORTSIZE)
    {
        for (i = 0; i < NSAMPLES; i++)
        {
            samplevals[i] = coord[arr[rand() % n]];
        }
        select_sort_dsample(samplevals, NSAMPLES);
        select_split(arr, n, samplevals[(NSAMPLES - 1) / 2], &start, &end,
                     coord);
        if (start > mid)
        {
            n = start;
        }
        else if (end <= mid)
        {
            arr += end;
            n -= end;
            mid -= end;
        }
        else
        {
            return;
        }
    }

    select_sort(arr, n, coord);
    return;
}

static void
select_split(int *arr, int n, double v, int *start, int *end, double *coord)
{
    int i, j, k;
    int t;

    i = 0;
    j = k = n;

    while (i < j)
    {
        if (coord[arr[i]] < v)
        {
            i++;
        }
        else if (coord[arr[i]] == v)
        {
            j--;
            SWAP(arr[i], arr[j], t);
        }
        else
        {
            j--;
            k--;
            t      = arr[i];
            arr[i] = arr[j];
            arr[j] = arr[k];
            arr[k] = t;
        }
    }
    *start = j;
    *end   = k;
    return;
}

static void
select_sort(int *arr, int n, double *coord)
{
    int i, j;
    int t;

    for (i = 1; i < n; i++)
    {
        t = arr[i];
        for (j = i; j > 0 && coord[arr[j - 1]] > coord[t]; j--)
        {
            arr[j] = arr[j - 1];
        }
        arr[j] = t;
    }
}

static void
select_sort_dsample(double *samp, int n)
{
    int i, j;
    double t;

    for (i = 1; i < n; i++)
    {
        t = samp[i];
        for (j = i; j > 0 && samp[j - 1] > t; j--)
        {
            samp[j] = samp[j - 1];
        }
        samp[j] = t;
    }
}
