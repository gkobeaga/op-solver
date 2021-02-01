#include "lp/lp.h"
#include "op-solver.h"

lp_data *
lp_create_data(int ncols, int nrows, int nzcnt)
{
    lp_data *data = malloc(sizeof(lp_data));

    data->name = NULL;

    if (nrows + ncols < 1)
        goto cleanup3;

    if (ncols > 0)
    {
        data->obj = malloc(ncols * sizeof(double));
        data->lb  = malloc(ncols * sizeof(double));
        data->ub  = malloc(ncols * sizeof(double));
        data->beg = malloc(ncols * sizeof(int));
        data->cnt = malloc(ncols * sizeof(int));

        if (!data->obj || !data->lb || !data->ub || !data->beg || !data->cnt)
        {
            fprintf(stderr, "%s\n", "lp_create_data: out of memory");
            goto cleanup1;
        }
    }
    else
    {
        data->obj = NULL;
        data->lb  = NULL;
        data->ub  = NULL;
    }

    if (nrows > 0)
    {
        data->rhs   = malloc(nrows * sizeof(double));
        data->sense = malloc(nrows * sizeof(char));
        if (!ncols)
        {
            data->beg = malloc(nrows * sizeof(int));
            data->cnt = malloc(nrows * sizeof(int));
        }
        if (!data->rhs || !data->sense)
        {
            fprintf(stderr, "%s\n", "lp_create_data: out of memory");
            goto cleanup2;
        }
    }
    else
    {
        data->rhs   = NULL;
        data->sense = NULL;
    }

    if (nzcnt > 0)
    {
        data->ind = malloc(nzcnt * sizeof(int));
        data->val = malloc(nzcnt * sizeof(double));
        if (!data->ind || !data->val)
        {
            fprintf(stderr, "%s\n", "lp_create_data: out of memory");
            goto cleanup3;
        }
    }
    else
    {
        data->ind = NULL;
        data->val = NULL;
    }

    data->ncols    = ncols;
    data->nrows    = nrows;
    data->nzcnt    = nzcnt;
    data->rowspace = nrows;
    data->colspace = ncols;
    data->nzspace  = nzcnt;

    return data;

cleanup3:
    if (nzcnt > 0)
    {
        free(data->ind);
        free(data->val);
    }
cleanup2:
    if (nrows > 0)
    {
        free(data->rhs);
        free(data->sense);
    }
cleanup1:
    if (ncols > 0)
    {
        free(data->obj);
        free(data->lb);
        free(data->ub);
    }

    data->ncols = 0;
    data->nrows = 0;
    data->nzcnt = 0;

    return NULL;
}

int
lp_realloc_data(lp_data *data, int count, double scale)
{
    int newsize = (int)(((double)data->nzspace) * scale);
    void *pind;
    void *pval;

    if (newsize < data->nzspace + 1000)
        newsize = data->nzspace + 1000;
    if (newsize < count)
        newsize = count;

    pind = realloc(data->ind, newsize * sizeof(int));
    pval = realloc(data->val, newsize * sizeof(double));
    if (!pind || !pval)
    {
        return 1;
    }
    else
    {
        data->ind     = pind;
        data->val     = pval;
        data->nzspace = newsize;
        return 0;
    }
}

void
lp_free_data(lp_data **data)
{
    if (*data)
    {
        if ((*data)->obj)
            free((*data)->obj);
        if ((*data)->ub)
            free((*data)->ub);
        if ((*data)->lb)
            free((*data)->lb);
        if ((*data)->rhs)
            free((*data)->rhs);
        if ((*data)->beg)
            free((*data)->beg);
        if ((*data)->cnt)
            free((*data)->cnt);
        if ((*data)->sense)
            free((*data)->sense);
        if ((*data)->ind)
            free((*data)->ind);
        if ((*data)->val)
            free((*data)->val);
        free(*data);
        *data = NULL;
    }
}
