#include "op-solver.h"

static void
__graph_create(solver_graph *graph),
__graph_delete(solver_graph *graph);

static void
__graph_create(solver_graph *graph)
{
    graph->nv_space  = 500;
    graph->na_space  = 100 * graph->nv_space;
    graph->nv        = 0;
    graph->n3v       = 0;
    graph->na        = 0;
    graph->v         = malloc(graph->nv_space * sizeof(graph_vertex *));
    graph->fixed     = NULL;
    graph->marker    = 0;
    graph->directed  = 0;
    graph->connected = 0;
    graph->archash   = NULL;
    graph->arcs      = malloc(graph->na_space * sizeof(graph_arc *));
    graph->data      = NULL;
    graph->shrunk    = NULL;
    graph->orig      = NULL;
    graph->tail      = NULL;
    graph->head      = NULL;
    graph_init_arc_hash(graph, 100);
    return;
}

solver_graph *
graph_create(void)
{
    solver_graph *graph = malloc(sizeof(solver_graph));
    __graph_create(graph);
    return graph;
}

int
graph_add_vertices(solver_graph *graph, int nadd)
{
    int i, nv_old;
    assert(nadd >= 0);
    assert(graph->nv + nadd < NV_MAX);

    nv_old = graph->nv;
    if (graph->nv_space < graph->nv + nadd)
    {
        realloc_scale(graph->v, graph->nv_space, graph->nv + nadd, 2);
    }

    for (i = nv_old; i < nv_old + nadd; i++)
    {
        graph_vertex *v;
        graph->v[i] = v = malloc(sizeof(graph_vertex));
        v->i            = i;
        v->mark         = 0;
        v->fixed        = 0;
        v->fixed_next   = NULL;
        v->fixed_prev   = NULL;
        v->branch       = 0;
        v->ind          = i;
        v->deg          = 0;
        v->coef         = 0.0;
        v->y            = 0.0;
        v->obj          = 0.0;
        v->vt_x         = 0.0;
        v->ut_x         = 0.0;
        v->comp         = 0;
        v->used         = 0;
        v->onecnt       = 0;
        v->onqueue      = 0;
        v->orig         = v;
        v->edge         = NULL;
        v->shrunk       = NULL;

        v->members  = NULL;
        v->nmembers = 1;
        v->parent   = v;
        v->max      = v;

        v->mark_aux = 0;

        v->depth = v->lowdepth = 0;
        v->active              = 0;

        v->prev = NULL;
        v->next = NULL;
        graph->nv++;
    }

    return nv_old;
}

graph_arc *
graph_add_arc(solver_graph *graph, int i, int j)
{
    graph_arc *arc, *rev;
    if (!(0 <= i && i < graph->nv))
    {
        printf("graph_add_arc: i = %d; tail vertex number out of range\n", i);
        exit(1);
    }
    if (!(0 <= j && j < graph->nv))
    {
        printf("graph_add_arc: j = %d; head vertex number out of range\n", j);
        exit(1);
    }

    if (!graph->directed)
    {
        if (i > j)
        {
            int temp;
            SWAP(i, j, temp);
        }
    }

    arc = graph_find_arc_hash(graph->archash, i, j);
    if (!arc)
    {
        if (graph->na == NA_MAX)
        {
            printf("graph_add_arc: too many arcs\n");
            exit(1);
        }
        arc       = malloc(sizeof(graph_arc));
        arc->tail = graph->v[i];
        arc->head = graph->v[j];

        arc->cost   = 0;
        arc->x      = 0.0;
        arc->obj    = 0.0;
        arc->age    = 0;
        arc->branch = 0;
        arc->fixed  = 0;
        arc->aux    = 0;

        arc->flow = 0.0;

        arc->used = 0;
        arc->coef = 0;

        if (!graph->v[i]->deg)
        {
            if (i == 0)
            {
                if (graph->tail)
                    graph->tail->prev = graph->v[i];
                else
                    graph->head = graph->v[i];
                graph->v[i]->next = graph->tail;
                graph->tail       = graph->v[i];
            }
            else
            {
                if (graph->head)
                    graph->head->next = graph->v[i];
                else
                    graph->tail = graph->v[i];
                graph->v[i]->prev = graph->head;
                graph->head       = graph->v[i];
            }
            graph->n3v++;
        }
        if (!graph->v[j]->deg)
        {
            if (j == 0)
            {
                graph->v[j]->next = graph->tail;
                graph->tail->prev = graph->v[j];
                graph->tail       = graph->v[j];
            }
            else
            {
                graph->v[j]->prev = graph->head;
                graph->head->next = graph->v[j];
                graph->head       = graph->v[j];
            }
            graph->n3v++;
        }

        graph_add_arc_hash(graph->archash, arc);

        if (graph->directed)
        {
            arc->t_prev = NULL;
            arc->t_next = graph->v[i]->out;
            if (arc->t_next != NULL)
                arc->t_next->t_prev = arc;

            arc->h_prev = NULL;
            arc->h_next = graph->v[j]->in;
            if (arc->h_next != NULL)
                arc->h_next->h_prev = arc;

            graph->v[i]->out = graph->v[j]->in = arc;
            rev = graph_find_arc_hash(graph->archash, j, i);
            if (rev)
            {
                arc->rev = rev;
                rev->rev = arc;
            }
        }
        else
        {
            arc->t_prev = NULL;
            arc->t_next = graph->v[i]->edge;
            if (arc->t_next != NULL)
            {
                if (arc->t_next->tail->i == arc->tail->i)
                    arc->t_next->t_prev = arc;
                else
                    arc->t_next->h_prev = arc;
            }

            arc->h_prev = NULL;
            arc->h_next = graph->v[j]->edge;
            if (arc->h_next != NULL)
            {
                if (arc->h_next->head->i == arc->head->i)
                    arc->h_next->h_prev = arc;
                else
                    arc->h_next->t_prev = arc;
            }

            graph->v[i]->edge = graph->v[j]->edge = arc;
        }

        if (graph->na + 1 >= graph->na_space)
        {
            realloc_scale(graph->arcs, graph->na_space, graph->na + 1, 1.3);
        }

        graph->arcs[graph->na] = arc;
        arc->ind               = graph->na;
        arc->i                 = graph->na;
        graph->na++;
        graph->v[i]->deg++;
        graph->v[j]->deg++;
    }

    return arc;
}

void
graph_del_arc(solver_graph *graph, graph_arc **arc)
{
    int rval = 0;
    graph_vertex *tail, *head;
    assert(graph->na > 0);
    assert(0 <= (*arc)->tail->i && (*arc)->tail->i < graph->nv);
    assert((*arc)->tail == graph->v[(*arc)->tail->i]);
    assert(0 <= (*arc)->head->i && (*arc)->head->i < graph->nv);
    assert((*arc)->head == graph->v[(*arc)->head->i]);
    tail = graph->v[(*arc)->tail->i];
    head = graph->v[(*arc)->head->i];

    if (graph_find_arc_hash(graph->archash, (*arc)->tail->i, (*arc)->head->i))
    {
        if ((*arc)->h_prev == NULL)
            (*arc)->head->edge = (*arc)->h_next;
        else
        {
            if ((*arc)->h_prev->head->i == (*arc)->head->i)
                (*arc)->h_prev->h_next = (*arc)->h_next;
            else
                (*arc)->h_prev->t_next = (*arc)->h_next;
        }
        if ((*arc)->h_next == NULL)
            ;
        else
        {
            if ((*arc)->h_next->head->i == (*arc)->head->i)
            {
                if ((*arc)->h_prev)
                    (*arc)->h_next->h_prev = (*arc)->h_prev;
                else
                    (*arc)->h_next->h_prev = NULL;
            }
            else
                (*arc)->h_next->t_prev = (*arc)->h_prev;
        }

        if ((*arc)->t_prev)
        {
            if ((*arc)->t_prev->tail->i == (*arc)->tail->i)
                (*arc)->t_prev->t_next = (*arc)->t_next;
            else
                (*arc)->t_prev->h_next = (*arc)->t_next;
        }
        else
            (*arc)->tail->edge = (*arc)->t_next;
        if ((*arc)->t_next)
        {
            if ((*arc)->t_next->tail->i == (*arc)->tail->i)
                (*arc)->t_next->t_prev = (*arc)->t_prev;
            else
                (*arc)->t_next->h_prev = (*arc)->t_prev;
        }
        else
        {
            if ((*arc)->t_prev != NULL)
            {
                if ((*arc)->t_prev->tail->i == (*arc)->tail->i)
                    (*arc)->t_prev->t_next = NULL;
                else
                    (*arc)->t_prev->h_next = NULL;
            }
        }

        rval = graph_del_arc_hash(graph->archash, (*arc));
        check_rval(rval, "arc hash not found\n", CLEANUP);

        graph->na--;

        int i;
        if ((*arc)->i < graph->na && graph->na > 0)
        {
            graph_arc *tmp;
            for (i = (*arc)->i; i < graph->na; i++)
            {
                tmp            = graph->arcs[i + 1];
                tmp->i         = i;
                tmp->ind       = i;
                graph->arcs[i] = tmp;
            }
        }

        graph->arcs[graph->na] = NULL;

        tail->deg--;
        head->deg--;
        if (!tail->deg)
        {
            if (tail->next)
                tail->next->prev = tail->prev;
            else
                graph->head = tail->prev;
            if (tail->prev)
                tail->prev->next = tail->next;
            else
                graph->tail = tail->next;
            graph->n3v--;
            tail->prev = NULL;
            tail->next = NULL;
        }
        if (!head->deg)
        {
            if (head->next)
                head->next->prev = head->prev;
            else
                graph->head = head->prev;
            if (head->prev)
                head->prev->next = head->next;
            else
                graph->tail = head->next;
            graph->n3v--;
            head->prev = NULL;
            head->next = NULL;
        }
        (*arc)->tail->y -= (*arc)->x / 2.0;
        (*arc)->head->y -= (*arc)->x / 2.0;

        free((*arc));
    }
    else
    {
        printf("not found %d %d\n", (*arc)->tail->i, (*arc)->head->i);
        exit(1);
    }

CLEANUP:
    return;
}

void
graph_delall_arcs(solver_graph *graph)
{
    int i;
    graph_arc *arc;

    if (!graph->na)
        return;

    // TODO: Optimize
    for (i = graph->na; i > 0; i--)
    {
        arc = graph->arcs[i - 1];
        graph_del_arc(graph, &arc);
    }

    free(graph->archash->table);
    free(graph->archash);
    graph_init_arc_hash(graph, 100);

    for (i = 0; i < graph->nv; i++) graph->v[i]->edge = NULL;

    free(graph->arcs);
    graph->arcs = malloc(graph->na_space * sizeof(graph_arc *));
    graph->na   = 0;
}

static void
__graph_delete(solver_graph *graph)
{
    int i;
    if (graph->archash)
        graph_free_arc_hash(&(graph->archash));
    for (i = 0; i < graph->nv; i++) free(graph->v[i]);
    free(graph->v);
    for (i = 0; i < graph->na; i++) free(graph->arcs[i]);
    free(graph->arcs);
    if (graph->shrunk)
    {
        graph_free(&(graph->shrunk));
    }
    return;
}

void
graph_erase(solver_graph *graph)
{
    __graph_delete(graph);
    __graph_create(graph);
    return;
}

void
graph_free(solver_graph **graph)
{
    if (*graph)
    {
        __graph_delete(*graph);
        free(*graph);
        *graph = NULL;
    }
    return;
}

graph_arc *
graph_find_arc(solver_graph *graph, graph_vertex *tail, graph_vertex *head)
{
    int end0, end1, t;
    end0 = tail->i;
    end1 = head->i;

    if (!graph->directed)
    {
        if (end0 > end1)
            SWAP(end0, end1, t);
    }

    return graph_find_arc_hash(graph->archash, end0, end1);
}

int
graph_copy(solver_graph *ingraph, solver_graph *outgraph)
{
    int rval = 0;
    int i;
    graph_vertex *v1, *v2;
    graph_arc *e;
    graph_arc *arc;

    outgraph->orig  = ingraph;
    ingraph->shrunk = outgraph;

    graph_add_vertices(outgraph, ingraph->nv);

    for (i = 0; i < ingraph->nv; i++)
    {
        outgraph->v[i]->fixed  = ingraph->v[i]->fixed;
        outgraph->v[i]->branch = ingraph->v[i]->branch;
        outgraph->v[i]->obj    = ingraph->v[i]->obj;
        outgraph->v[i]->orig   = ingraph->v[i];
        ingraph->v[i]->shrunk  = outgraph->v[i];
    }

    outgraph->fixed = NULL;
    for (v1 = ingraph->fixed; v1; v1 = v1->fixed_next)
    {
        if (outgraph->fixed)
            outgraph->fixed->fixed_prev = outgraph->v[v1->i];
        outgraph->v[v1->i]->fixed_next = outgraph->fixed;
        outgraph->fixed                = outgraph->v[v1->i];
    }

    (outgraph->marker)++;
    outgraph->n3v = 0;
    for (i = 0; i < ingraph->na; i++)
    {
        arc = ingraph->arcs[i];
        {
            e  = graph_add_arc(outgraph, arc->tail->i, arc->head->i);
            v1 = outgraph->v[arc->tail->i];
            v2 = outgraph->v[arc->head->i];

            e->x    = arc->x;
            e->cost = arc->cost;
            e->obj  = arc->obj;
            e->tail->y += e->x / 2.0;
            e->head->y += e->x / 2.0;

            if (arc->x > SOLVER_ONEMINUS)
            {
                v1->onecnt++;
                v2->onecnt++;
                outgraph->onecnt++;
            }
        }
    }
    return rval;
}

void
graph_print(solver_graph *graph)
{
    int i;
    graph_vertex *v, *other;
    graph_arc *e, *enext;

    printf("\n### Support Graph\n");

    for (i = 0; i < graph->nv; i++)
    {
        v = graph->v[i];
        if (v->y)
        {
            printf("    Node %d [%.3f]:", v->i, v->y);
            for (e = v->edge; e; e = enext)
            {
                enext = outnext(e, v);
                other = otherend(e, v);
                printf("%d [%.3f] ", other->i, e->x);
            }
            printf("\n");
        }
    }

    printf("### End Graph\n\n");
}

void
graph_write(solver_graph *graph, char *fname)
{
    int i;
    FILE *file;
    graph_vertex *v, *other;
    graph_arc *e, *enext;

    printf("\n");
    printf("Writing graph to '%s'...\n", fname);

    file = fopen(fname, "w");
    if (file == NULL)
    {
        printf("Unable to create '%s'\n", "saved.graph");
        exit(1);
    }

    fprintf(file, "Number_of_nodes: %d\n", graph->nv);
    fprintf(file, "Number_of_arcs: %d\n", graph->na);
    for (i = 0; i < graph->nv; i++)
    {
        v = graph->v[i];
        if (v->y)
        {
            fprintf(file, "Node %d [%.6f]:", v->i, v->y);
            for (e = v->edge; e; e = enext)
            {
                enext = outnext(e, v);
                other = otherend(e, v);
                if (e->tail->i == v->i)
                    fprintf(file, " %d [%.6f]", other->i,
                            roundf(e->x * ((double)10000000)) /
                            ((double)10000000));
            }
            fprintf(file, "\n");
        }
    }

    if (file != NULL)
        fclose(file);
}

solver_graph *
graph_read(char *fname)
{
    int rval   = 0;
    FILE *file = NULL;
    int nv, na;
    graph_arc *arc;
    char buf[256];
    solver_graph *graph = graph_create();

    file = fopen(fname, "r");
    if (file == NULL)
    {
        printf("Unable to open '%s'\n", fname);
        goto done;
    }

    rval = fscanf(file, "%*s %d", &nv);
    check_assert(rval == 1, "failed", done);
    rval = fscanf(file, "%*s %d", &na);
    check_assert(rval == 1, "failed", done);

    graph_add_vertices(graph, nv);

    int tail, head, n_read = 0;
    double node_val, arc_val;
    char *read_ptr = buf;
    while (fgets(buf, 254, file) != NULL)
    {
        if (sscanf(buf, "%*s %d [%lf]: %n", &tail, &node_val, &n_read) > 0)
        {
            graph->v[tail]->y = node_val;
            read_ptr += n_read;
        }
        while (sscanf(read_ptr, "%d [%lf] %n", &head, &arc_val, &n_read) == 2)
        {
            arc    = graph_add_arc(graph, tail, head);
            arc->x = arc_val;
            read_ptr += n_read;
        }
        read_ptr = buf;
    }

    if (file != NULL)
        fclose(file);

    return graph;

done:
    return NULL;
}

void
graph_identify_vertices(solver_graph *graph, graph_vertex *v, graph_vertex *u)
{
    double oldobj, oldcost;
    graph_arc *e, *new, *old, *next;
    graph_vertex *other;

    assert(v->i != u->i);

    u->parent = v;

    if (!v->members)
    {
        v->members = u;
    }
    else if (!u->members)
    {
        u->members = v->members;
        v->members = u;
    }
    else
    {
        graph_vertex *t;
        for (t = v->members; t->members; t = t->members)
            ;
        t->members = u;
    }
    v->nmembers += u->nmembers;
    v->obj += u->obj;

    if (u->fixed)
    {
        if (!v->fixed)
        {
            if (u->fixed_prev)
                u->fixed_prev->fixed_next = v;
            else
                graph->fixed = v;
            if (u->fixed_next)
                u->fixed_next->fixed_prev = v;
            v->fixed_next = u->fixed_next;
            v->fixed_prev = u->fixed_prev;
            v->branch     = u->branch;
        }
        else
        {
            if (u->fixed_prev)
                u->fixed_prev->fixed_next = u->fixed_next;
            else
                graph->fixed = u->fixed_next;
            if (u->fixed_next)
                u->fixed_next->fixed_prev = u->fixed_prev;
        }
        v->fixed      = u->fixed;
        u->fixed_next = NULL;
        u->fixed_prev = NULL;
        if (u->branch < v->branch)
            v->branch = u->branch;
    }

    if (v->max->orig->y < u->max->orig->y)
        v->max = u->max;

    old     = graph_find_arc_hash(graph->archash, v->i, u->i);
    oldobj  = old ? old->obj : 0;
    oldcost = old ? old->cost : 0;
    for (e = u->edge; e; e = next)
    {
        other = otherend(e, u);
        next  = outnext(e, u);
        if (other->i != v->i)
        {
            new = graph_add_arc(graph, v->i, other->i);
            new->flow += e->flow;

            if (new->x > SOLVER_ONEMINUS)
            {
                new->tail->onecnt--;
                new->head->onecnt--;
                graph->onecnt--;
            }
            new->x += e->x;
            new->obj  = e->obj + u->obj + oldobj;
            new->cost = e->cost + oldcost;
            v->y += e->x / 2.0;
            other->y += e->x / 2.0;
            if (new->x > SOLVER_ONEMINUS)
            {
                new->tail->onecnt++;
                new->head->onecnt++;
                graph->onecnt++;
            }
        }

        if (e->x > SOLVER_ONEMINUS)
        {
            e->tail->onecnt--;
            e->head->onecnt--;
            graph->onecnt--;
        }
        graph_del_arc(graph, &e);
    }
}

int
graph_expand_vertex(solver_graph *graph, graph_vertex *srkv, int *vcount,
                    graph_vertex ***verts)
{
    int rval = 0;
    int cnt;
    graph_vertex **tverts;
    graph_vertex *v;

    *vcount = 0;
    *verts  = NULL;

    cnt = 1;
    for (v = srkv->members; v; v = v->members) cnt++;

    tverts = malloc(cnt * sizeof(graph_vertex *));
    check_null(tverts, "graph_expand_vertex failed", CLEANUP);

    tverts[0] = srkv->orig;
    cnt       = 1;
    for (v = srkv->members; v; v = v->members) tverts[cnt++] = v->orig;

    *vcount = cnt;
    *verts  = tverts;

CLEANUP:

    return rval;
}

static int
sort_nodes(const void *vv, const void *uu);
int
graph_reorder_vertices(solver_graph *graph)
{
    int i, k;

    if (!graph->nv)
        return 0;

    graph_vertex *v;

    graph_vertex **nnvertices = NULL;

    // Sort nodes
    nnvertices = malloc(graph->n3v * sizeof(graph_vertex *));
    for (v = graph->tail, k = 0; v; v = v->next) nnvertices[k++] = v;
    assert(k == graph->n3v);
    qsort(nnvertices, k, sizeof(graph_vertex *), sort_nodes);

    graph->tail = nnvertices[0];
    v           = nnvertices[0];
    v->prev     = NULL;
    v->mark_aux = 1;
    for (i = 1; i < k; i++)
    {
        v->next             = nnvertices[i];
        nnvertices[i]->prev = v;
        v                   = nnvertices[i];
        v->next             = NULL;
        v->mark_aux         = 1;
    }
    graph->head = v;

    if (nnvertices)
        free(nnvertices);

    return 0;
}

static int
sort_nodes(const void *vv, const void *uu)
{
    graph_vertex *v = *(graph_vertex **)vv, *u = *(graph_vertex **)uu;

#if REORDER == 0
    if (v->fixed && u->fixed)
    {
        if (rand() % 2)
            return -1;
        else
            return +1;
    }
#elif REORDER == 1
    if (v->fixed && u->fixed)
    {
        if (v->i < u->i)
            return -1;
        if (v->i > u->i)
            return +1;
    }
#else
    if (v->fixed && u->fixed)
        return 0;
#endif
    if (v->fixed)
        return -1;
    if (u->fixed)
        return +1;

    if (v->y + SOLVER_ZEROPLUS < u->y)
        return +1;
    if (v->y - SOLVER_ZEROPLUS > u->y)
        return -1;

#if REORDER == 0
    if (rand() % 2)
        return -1;
    else
        return +1;
#elif REORDER == 1
    if (v->i < u->i)
        return -1;
    if (v->i > u->i)
        return +1;
#else
    return 0;
#endif
}

void
graph_plot(solver_graph *graph)
{
    int i;
    graph_arc *arc;

    double *x = graph->data->map->x;
    double *y = graph->data->map->y;

    FILE *pipe = popen("gnuplot -persist", "w");

    fprintf(pipe, "set terminal pdfcairo\n");
    fprintf(pipe, "set output 'graph.pdf'\n");
    fprintf(pipe, "set multiplot\n");
    fprintf(pipe, "unset ytics\n");
    fprintf(pipe, "unset xtics\n");
    fprintf(pipe, "set title 'CP Solution'\n");
    fprintf(pipe, "set xlabel 'X'\n");
    fprintf(pipe, "set ylabel 'Y'\n");
    fprintf(pipe,
            "set offset graph 0.05, graph 0.05, graph 0.05, graph 0.05\n");

    fprintf(pipe, "set nokey\n");

    // Arcs
    fprintf(pipe, "set style arrow 1 nohead lt 1 lc rgb 'black'\n");
    fprintf(pipe, "set style arrow 2 nohead lt 1 lc rgb 'red'\n");
    fprintf(pipe, "set style arrow 3 nohead lt 0 lw 2\n");
    for (i = 0; i < graph->na; i++)
    {
        arc = graph->arcs[i];
        if (arc->x > SOLVER_ONEMINUS)
            fprintf(pipe, "set arrow from %lf,%lf to %lf,%lf as 1\n",
                    x[arc->tail->i], y[arc->tail->i], x[arc->head->i],
                    y[arc->head->i]);
        else if (arc->x > SOLVER_ZEROPLUS && arc->x > 0.5)
            fprintf(pipe, "set arrow from %lf,%lf to %lf,%lf as 2\n",
                    x[arc->tail->i], y[arc->tail->i], x[arc->head->i],
                    y[arc->head->i]);
        else if (arc->x > SOLVER_ZEROPLUS && arc->x < 0.5)
            fprintf(pipe, "set arrow from %lf,%lf to %lf,%lf as 3\n",
                    x[arc->tail->i], y[arc->tail->i], x[arc->head->i],
                    y[arc->head->i]);
        else if (arc->x > SOLVER_ZEROPLUS && arc->x == 0.5)
            fprintf(pipe, "set arrow from %lf,%lf to %lf,%lf as 2\n",
                    x[arc->tail->i], y[arc->tail->i], x[arc->head->i],
                    y[arc->head->i]);
    }

    // Points
    fprintf(pipe, "set style fill solid 1.0 border -1\n");
    fprintf(pipe, "set style circle radius screen 0.007\n");

    fprintf(pipe, "plot '-' w circles fc rgb 'green', ");
    fprintf(pipe, "     '-' w circles fc rgb 'black', ");
    fprintf(pipe, "     '-' w circles fc rgb 'red', ");
    fprintf(pipe, "     '-' w circles fc rgb 'gray', ");
    fprintf(pipe, "     '-' w labels tc rgb 'white' font 'Times Roman,5', ");
    fprintf(pipe, "     '-' w labels tc rgb 'black' font 'Times Roman,5'\n");

    fprintf(pipe, "%lf %lf\n", x[0], y[0]);
    fprintf(pipe, "e\n");

    for (i = 1; i < graph->nv; i++)
    {
        if (graph->v[i]->y > SOLVER_ONEMINUS)
            fprintf(pipe, "%lf %lf as 1\n", x[i], y[i]);
    }
    fprintf(pipe, "e\n");

    for (i = 1; i < graph->nv; i++)
    {
        if (graph->v[i]->y < SOLVER_ONEMINUS && graph->v[i]->y >= 0.5)
            fprintf(pipe, "%lf %lf as 1\n", x[i], y[i]);
    }
    fprintf(pipe, "e\n");

    for (i = 1; i < graph->nv; i++)
    {
        if (graph->v[i]->y < 0.5 && graph->v[i]->y > SOLVER_ZEROPLUS)
            fprintf(pipe, "%lf %lf as 1\n", x[i], y[i]);
    }
    fprintf(pipe, "e\n");

    for (i = 0; i < graph->nv; i++)
    {
        if (i != 0 && graph->v[i]->y > SOLVER_ONEMINUS)
            fprintf(pipe, "%lf %lf %d\n", x[i], y[i], i);
    }
    fprintf(pipe, "e\n");

    for (i = 0; i < graph->nv; i++)
    {
        if (i == 0 || graph->v[i]->y <= SOLVER_ONEMINUS)
            fprintf(pipe, "%lf %lf %d\n", x[i], y[i], i);
    }
    fprintf(pipe, "e\n");

    fclose(pipe);

    printf("Press Enter to Continue\n");
    while (getchar() != '\n')
        ;

    return;
}
