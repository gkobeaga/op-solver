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
    FILE *file = NULL;

    if (sol == NULL)
    {
        printf("No OP solution to write.\n");
    }

    printf("\n");
    printf("Writing OP solution to '%s'...\n", fname);

    file = fopen(fname, "a");
    check_null(file, "Unable to create file", done);

    fprintf(file, "{ ");

    // Problem
    fprintf(file, "\"prob\": { ");
    fprintf(file, "\"name\": \"%s\", ", cp->data->name);
    fprintf(file, "\"type\": \"OP\",");
    fprintf(file, "\"n\": %d, ", cp->n);
    fprintf(file, "\"d0\": %.0f, ", cp->data->cap);
    fprintf(file, "\"depot\": %d", cp->from + 1);
    fprintf(file, "}, ");

    // Solution
    fprintf(file, "\"sol\": { ");
    fprintf(file, "\"val\": %.0f, ", cp->sol->val);
    fprintf(file, "\"cap\": %.0f, ", cp->sol->cap);
    fprintf(file, "\"sol_ns\": %d, ", cp->sol->ns);
    if (cp->ip && cp->ip->sol)
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
    for (i = 0; i < sol->ns - 1; i++) fprintf(file, "%d, ", sol->cycle[i] + 1);
    fprintf(file, "%d]", sol->cycle[sol->ns - 1] + 1);
    fprintf(file, "}");

    fprintf(file, "}");

    rval = 0;
done:

    if (file != NULL)
        fclose(file);

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

cp_sol *
cp_get_sol_from_cycle(cp_prob *cp, int ns, int *cycle)
{
    double cap = 0, val = 0;

    cp_sol *sol = cp_create_sol(cp);

    sol->selected[cycle[0]] = 1;
    sol->sposition[0]       = cycle[0];
    sol->cycle[0]           = cycle[0];
    cap += data_get_norm(cp->data, cycle[ns - 1], cycle[0]);
    val += cp->data->obj_node[cycle[0]];

    for (int i = 1; i < ns; i++)
    {
        int prev = cycle[i - 1];
        int cur  = cycle[i];

        sol->cod_fr[prev]  = cur;
        sol->cod_bk[cur]   = prev;
        sol->selected[cur] = 1;
        sol->sposition[i]  = cur;
        sol->cycle[i]      = cur;

        cap += data_get_norm(cp->data, prev, cur);
        val += cp->data->obj_node[cur];
    }

    for (int i = 0, j = 0; i < cp->n; i++)
    {
        if (sol->selected[i])
            sol->sposition[j++] = i;
    }

    sol->val = val;
    sol->cap = cap;
    sol->ns  = ns;
    return sol;
}

void
cp_plot_sol(cp_prob *cp, cp_sol *sol)
{
    solver_graph *graph = cp_conv_sol_to_graph(cp, sol);

    graph_plot(graph);

    graph_free(&graph);
}
