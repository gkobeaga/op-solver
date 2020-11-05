#include "op-solver.h"

/*  This implementation of the Push-Relabel Flow Algorithm described in  */
/*  A. Goldberg and R. Tarjan, "A new approach to the maximum-flow       */
/*  problem"  is based on Concorde.                                      */

static void
set_labels(solver_graph *graph, graph_vertex *s, graph_vertex *t,
           graph_vertex **high, int *highest, graph_vertex **level),
relabel_dist(int n3v, graph_vertex **level, graph_vertex *v),
add_to_active(graph_vertex **high, int *highest, graph_vertex *v),
backwards_bfs(solver_graph *graph, graph_vertex *t, graph_vertex **high,
              int *highest, graph_vertex **level, int marker);

double
graph_get_maxflow_st(solver_graph *graph, graph_vertex *s, graph_vertex *t)
{
    graph_arc *e;
    int count;
    int i;
    double residual_flow;
    graph_vertex **high;
    graph_vertex **level;
    int highest;
    graph_vertex *v;
    graph_vertex *other;

    level = malloc((graph->n3v + 1) * sizeof(graph_vertex *));
    for (i = 0; i < graph->n3v; i++) level[i] = NULL;
    level[graph->n3v] = NULL;

    high = malloc(graph->n3v * sizeof(graph_vertex *));

    i = 0;
    for (v = graph->tail; v; v = v->next)
    {
        v->excess = 0.0;
        v->active = 0;
        high[i++] = NULL;
    }
    highest = 0;

    for (i = 0; i < graph->na; i++) graph->arcs[i]->flow = 0.0;

    t->active = 1;
    s->active = 1;

    for (e = s->edge; e; e = outnext(e, s))
    {
        other = otherend(e, s);
        if (e->x > 0.0)
        {
            if (e->tail->i == s->i)
            {
                e->flow = e->x;
                other->excess += e->x;
            }
            else
            {
                e->flow = -e->x;
                other->excess += e->x;
            }
            other->active = 1;
        }
    }

    set_labels(graph, s, t, high, &highest, level);

    count = 0;
    while (highest)
    {
        v             = high[highest];
        v->active     = 0;
        high[highest] = high[highest]->high_next;
        if (!high[highest])
        {
            highest--;
            while (highest && (high[highest] == NULL)) highest--;
        }

        if (count == graph->n3v)
        {
            set_labels(graph, s, t, high, &highest, level);
            // If v is not reachable from sink, continue.
            if (v->label >= graph->n3v)
                continue;
            count = 0;
        }
        else
            count++;

        if (v->ecurrent)
        {
            while (v->excess > 0.0)
            {
                e     = v->ecurrent;
                other = otherend(e, v);

                if (e->tail->i == v->i)
                    residual_flow = e->x - e->flow;
                else
                    residual_flow = e->x + e->flow;

                if (v->label == other->label + 1 && residual_flow > 0.0)
                {
                    if (v->excess <= residual_flow)
                    {
                        if (e->tail->i == v->i)
                            e->flow += v->excess;
                        else
                            e->flow -= v->excess;
                        other->excess += v->excess;
                        v->excess = 0.0;

                        add_to_active(high, &highest, other);
                    }
                    else
                    {
                        if (e->tail->i == v->i)
                            e->flow += residual_flow;
                        else
                            e->flow -= residual_flow;
                        other->excess += residual_flow;
                        v->excess -= residual_flow;

                        add_to_active(high, &highest, other);
                        v->ecurrent = outnext(e, v);
                        if (!v->ecurrent)
                        {
                            v->ecurrent = v->edge;
                            relabel_dist(graph->n3v, level, v);
                            break;
                        }
                    }
                }
                else
                {
                    v->ecurrent = outnext(e, v);
                    if (!v->ecurrent)
                    {
                        v->ecurrent = v->edge;
                        relabel_dist(graph->n3v, level, v);
                        break;
                    }
                }
            }
        }
        else
        {
            relabel_dist(graph->n3v, level, v);
        }

        if (v->excess > 0.0 && v->label < graph->n3v)
        {
            add_to_active(high, &highest, v);
        }
    }

    free(high);
    free(level);

    return t->excess;
}

static void
add_to_active(graph_vertex **high, int *highest, graph_vertex *v)
{
    if (!v->active)
    {
        v->high_next   = high[v->label];
        high[v->label] = v;
        if (*highest < v->label)
            *highest = v->label;
        v->active = 1;
    }
}

static void
relabel_dist(int n3v, graph_vertex **level, graph_vertex *v)
{
    int newlabel = SOLVER_MAXINT;
    graph_arc *rele;
    graph_vertex *other;
    int tmplabel, oldlabel;

    for (rele = v->edge; rele; rele = outnext(rele, v))
    {
        other    = otherend(rele, v);
        tmplabel = other->label;
        if (rele->tail->i == v->i)
        {
            if (rele->x - rele->flow > 0.0 && tmplabel < newlabel)
                newlabel = tmplabel;
        }
        else
        {
            if (rele->x + rele->flow > 0.0 && tmplabel < newlabel)
                newlabel = tmplabel;
        }
    }

    oldlabel = v->label;
    v->label = ++newlabel;

    if (v->level_prev)
        v->level_prev->level_next = v->level_next;
    else
        level[oldlabel] = v->level_next;
    if (v->level_next)
        v->level_next->level_prev = v->level_prev;

    // If v is reachable from sink
    if (newlabel < n3v)
    { // Insert node in the new label level
        if (level[newlabel])
        {
            level[newlabel]->level_prev = v;
            v->level_next               = level[newlabel];
            v->level_prev               = NULL;
            level[newlabel]             = v;
        }
        else
        {
            v->level_prev   = NULL;
            v->level_next   = NULL;
            level[newlabel] = v;
        }

        // Drop node from the old label level
        if (!level[oldlabel])
        {
            oldlabel++;
            while (level[oldlabel])
            {
                graph_vertex *relno;
                for (relno = level[oldlabel]; relno; relno = relno->level_next)
                    relno->label = n3v;
                level[oldlabel] = NULL;
                oldlabel++;
            }
        }
    }
}

static void
set_labels(solver_graph *graph, graph_vertex *s, graph_vertex *t,
           graph_vertex **high, int *highest, graph_vertex **level)
{
    graph_vertex *v;
    int marker = ++(graph->marker);

    t->label = 0;
    backwards_bfs(graph, t, high, highest, level, marker);
    if (s->mark == marker)
    {
        printf("Help - s should not get a label\n");
        s->label = graph->n3v;
    }

    for (v = graph->tail; v; v = v->next)
    {
        v->ecurrent = v->edge;
        if (v->mark != marker)
            v->label = graph->n3v;
    }
}

// Backward breadth-first search from sink
static void
backwards_bfs(solver_graph *graph, graph_vertex *t, graph_vertex **high,
              int *highest, graph_vertex **level, int marker)
{
    int i;
    graph_vertex *v, *next, *other;
    graph_arc *e;
    int label;
    graph_vertex *dummy = malloc(sizeof(graph_vertex));

    t->mark = marker;
    label   = t->label;

    // Initialize: label and high vectors
    for (i = 0; level[i]; i++) level[i] = NULL;
    level[label]  = t;
    t->level_next = NULL;

    for (i = 0; i <= *highest; i++) high[i] = NULL;
    *highest = 0;

    // Start search from sink
    next           = t;
    t->search_next = NULL;

    // Perform search
    do
    { // Initialize new level
        level[label]->level_prev = NULL;
        level[label + 1]         = dummy;
        dummy->level_prev        = NULL;
        dummy->level_next        = NULL;

        label++;
        for (v = next, next = NULL; v; v = v->search_next)
        {
            for (e = v->edge; e; e = outnext(e, v))
            {
                other = otherend(e, v);
                if (e->tail->i == v->i)
                {
                    if (other->mark != marker && e->x + e->flow > 0.0)
                    {
                        other->label             = label;
                        other->search_next       = next;
                        next                     = other;
                        other->mark              = marker;
                        other->level_next        = level[label];
                        level[label]->level_prev = other;
                        level[label]             = other;
                        if (other->active)
                        {
                            other->high_next = high[label];
                            high[label]      = other;
                        }
                    }
                }
                else
                {
                    if (other->mark != marker && e->x - e->flow > 0.0)
                    {
                        other->label             = label;
                        other->search_next       = next;
                        next                     = other;
                        other->mark              = marker;
                        other->level_next        = level[label];
                        level[label]->level_prev = other;
                        level[label]             = other;
                        if (other->active)
                        {
                            other->high_next = high[label];
                            high[label]      = other;
                        }
                    }
                }
            }
        }
        if (dummy->level_prev)
        {
            dummy->level_prev->level_next = NULL;
            level[label]->level_prev      = NULL;
        }
        else
        {
            level[label] = NULL;
        }
        if (high[label])
            *highest = label;
    } while (next);

    free(dummy);
}
