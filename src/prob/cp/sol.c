#include "op-solver.h"

static void
cp_create_sol_(cp_sol *sol, int n)
{
    sol->tot_n     = n;
    sol->selected  = malloc(n * sizeof(int));
    sol->sposition = malloc(n * sizeof(int));
    sol->cycle     = malloc(n * sizeof(int));
    sol->cod_fr    = malloc(n * sizeof(int));
    sol->cod_bk    = malloc(n * sizeof(int));
    sol->val       = 0.0;
    sol->cap       = 0.0;
    sol->ns        = 0;
    memset(sol->selected, 0, n * sizeof(int));
    memset(sol->sposition, -1, n * sizeof(int));
    memset(sol->cycle, -1, n * sizeof(int));
    for (int i = 0; i < n; i++)
    {
        sol->cod_bk[i] = i;
        sol->cod_fr[i] = i;
    }
}

cp_sol *
cp_create_sol(cp_prob *cp)
{
    cp_sol *sol = malloc(sizeof(cp_sol));
    cp_create_sol_(sol, cp->n);

    return sol;
}

static void
cp_delete_sol_(cp_sol *sol)
{
    free(sol->selected);
    free(sol->sposition);
    free(sol->cycle);
    free(sol->cod_bk);
    free(sol->cod_fr);
}

void
cp_erase_sol(cp_sol *sol)
{
    int n = sol->tot_n;
    cp_delete_sol_(sol);
    cp_create_sol_(sol, n);
    return;
}

void
cp_free_sol(cp_sol **sol)
{
    cp_delete_sol_(*sol);
    free(*sol);
    *sol = NULL;
    return;
}

void
cp_copy_sol(cp_sol *insol, cp_sol *outsol)
{
    int i;
    cp_erase_sol(outsol);
    for (i = 0; i < insol->tot_n; i++)
    {
        if (insol->cod_fr)
            outsol->cod_fr[i] = insol->cod_fr[i];
        if (insol->cod_bk)
            outsol->cod_bk[i] = insol->cod_bk[i];
        if (insol->selected)
            outsol->selected[i] = insol->selected[i];
        if (insol->sposition)
            outsol->sposition[i] = insol->sposition[i];
        if (insol->cycle)
            outsol->cycle[i] = insol->cycle[i];
    }
    outsol->val = insol->val;
    outsol->cap = insol->cap;
    outsol->ns  = insol->ns;
    return;
}

int
cp_get_sol_from_graph(cp_prob *cp, solver_graph *graph, cp_sol **sol)
{
    int rval = 0;
    int i;
    graph_arc *arc, *next;
    graph_vertex *prev, *other;

    if (!(*sol))
    {
        *sol = cp_create_sol(cp);
        check_null(*sol, " out of memory", CLEANUP);
    }
    else
        cp_erase_sol(*sol);

    prev  = graph->v[0];
    other = NULL;
    for (arc = graph->v[0]->edge; arc && (!other || other->i); arc = next)
    {
        other = otherend(arc, prev);
        (*sol)->val += cp->data->obj_node[prev->i];
        (*sol)->cycle[(*sol)->ns++] = prev->i;
        (*sol)->cod_fr[prev->i]     = other->i;
        (*sol)->cod_bk[other->i]    = prev->i;
        (*sol)->selected[prev->i]   = 1;
        (*sol)->cap += data_get_norm(cp->data, prev->i, other->i);
        next = outnext(arc, other) ? outnext(arc, other) : other->edge;
        prev = other;
    }

    for (i = 0; i < graph->nv; i++)
    {
        (*sol)->sposition[i] = cp->n;
    }

    (*sol)->ns = 0;
    for (i = 0; i < graph->nv; i++)
    {
        if ((*sol)->selected[i])
        {
            (*sol)->sposition[(*sol)->ns] = i;
            (*sol)->ns++;
        }
    }

CLEANUP:

    return rval;
}

void
cp_print_sol(cp_prob *cp, cp_sol *sol)
{
    printf("CP solution:\n");
    printf(" - Total nodes: %d\n", sol->tot_n);
    printf(" - Visited nodes: %d\n", sol->ns);
    printf(" - Objetive value: %f\n", sol->val);
    printf(" - Cycle length: %f\n", sol->cap);
#if 0
    printf(" - Cycle: ");
    for (int i = 0; i < 3; i++)
        printf("%d ", cp->data->map->inv[sol->cycle[i]]);
    printf(" ... ");
    for (int i = 0; i < 3; i++)
        printf("%d ", cp->data->map->inv[sol->cycle[sol->ns - 3 + i]]);
    printf("\n");
#else
    printf(" - Cycle: ");
    for (int i = 0; i < sol->ns; i++)
        printf("%d ", cp->data->map->inv[sol->cycle[i]]);
    printf("\n");
#endif
}

int
cp_write_sol(cp_prob *cp, cp_sol *sol, const char *fname)
{
    int rval, i;
    FILE *fp = NULL;

    if (sol == NULL)
    {
        printf("No OP solution to write.\n");
    }

    printf("\n");
    printf("Writing OP solution to '%s'...\n", fname);

    fp = fopen(fname, "w");
    if (fp == NULL)
    {
        printf("Unable to create '%s'\n", fname);
        rval = 1;
        goto done;
    }

    fprintf(fp, "NAME : %s\n", cp->data->name);
    fprintf(fp, "TYPE : OP\n");
    fprintf(fp, "DIMENSION : %d\n", cp->n);
    fprintf(fp, "COST_LIMIT : %.2f\n", cp->cap);
    fprintf(fp, "ROUTE_NODES : %d\n", sol->ns);
    fprintf(fp, "ROUTE_SCORE : %.2f\n", sol->val);
    fprintf(fp, "ROUTE_COST : %.2f\n", sol->cap);
    fprintf(fp, "NODE_SEQUENCE_SECTION\n");
    for (i = 0; i < sol->ns; i++) fprintf(fp, "%d\n", sol->cycle[i] + 1);
    fprintf(fp, "-1\n");
    fprintf(fp, "DEPOT_SECTION\n");
    fprintf(fp, "%d\n", cp->from + 1);
    fprintf(fp, "-1\n");
    fprintf(fp, "EOF\n");

    rval = 0;
done:

    if (fp != NULL)
        fclose(fp);

    return rval;
}

solver_graph *
cp_conv_sol_to_graph(cp_prob *cp, cp_sol *sol)
{
    solver_graph *graph = graph_create();
    graph_arc *arc;
    graph_add_vertices(graph, sol->tot_n);
    graph->data = cp->data;

    for (int i = 0; i < sol->ns; i++)
    {
        arc = graph_add_arc(graph, sol->cycle[i], sol->cod_fr[sol->cycle[i]]);
        arc->x       = 1;
        arc->tail->y = 1;
        arc->head->y = 1;
    }

    return graph;
}

void
cp_plot_sol(cp_prob *cp, cp_sol *sol)
{
    solver_graph *graph = cp_conv_sol_to_graph(cp, sol);

    graph_plot(graph);

    graph_free(&graph);
}
