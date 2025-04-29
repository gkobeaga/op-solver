#include "op-solver.h"

/* Implementation of the Gomory-hu tree construction based on Concorde  */

static void
unshrink(solver_graph *graph, graph_vertex *srknode, int ecount, int *elist,
         double *ecap),
ghtree_free_node(graph_ghtree_node *p), update_outmax(graph_ghtree_node *n),
ghtree_print_work(graph_ghtree_node *root);

static int
shrinkdown(graph_ghtree *ghtree, solver_graph *graph, int nodecount,
           graph_vertex **nodeset, graph_vertex *srknode, int *ecount,
           int *elist, double *ecap),
gh_recursive(solver_graph *graph, graph_ghtree *ghtree, graph_ghtree_node *n,
             int scount, graph_vertex **special),
countdescendants(graph_ghtree_node *n);

int
graph_get_ghtree(solver_graph *graph, int rootind, graph_ghtree **ghtree)
{
    int rval = 0;
    graph_vertex **special;
    int scount;
    graph_vertex *v;
    double maxval = -SOLVER_MAXDOUBLE;
    graph_vertex *maxv;
    int tail_max = 1;

    for (v = graph->tail; v; v = v->next)
    {
        if (graph->tail->y < v->y - SOLVER_ZEROPLUS)
        {
            tail_max = 0;
        }
        if (v->max->y > maxval)
        {
            maxval = v->max->y;
            maxv   = v->max;
        }
        v->mark_aux = 1;
    }

    *ghtree = graph_create_ghtree();
    graph_add_ghtree_nodes(*ghtree, graph->nv);
    if (graph->v[rootind]->deg)
        (*ghtree)->root = (*ghtree)->n[rootind];
    else
    {
        printf("Not a valid Gomory-Hu tree root\n");
        exit(1);
    }
    (*ghtree)->root->cutval = SOLVER_MAXDOUBLE;

    special = malloc(graph->n3v * sizeof(graph_vertex *));
    scount  = 0;
    for (v = graph->tail; v; v = v->next)
    {
        if (v->mark_aux == 1)
            special[scount++] = v;
        (*ghtree)->n[v->i]->in_max  = v;
        (*ghtree)->n[v->i]->out_max = graph->tail;
        (*ghtree)->n[v->i]->special = v;
    }
    (*ghtree)->root->in_max  = maxv;
    (*ghtree)->root->out_max = special[0];
    (*ghtree)->root->special = special[0];

    rval = gh_recursive(graph, *ghtree, (*ghtree)->root, scount, special);
    check_rval(rval, "gh_recursive failed", CLEANUP);

    countdescendants((*ghtree)->root);
    if (!tail_max)
        update_outmax((*ghtree)->root);

CLEANUP:

    if (special)
        free(special);

    return rval;
}

static int
gh_recursive(solver_graph *graph, graph_ghtree *ghtree, graph_ghtree_node *n,
             int scount, graph_vertex **special)
{
    int rval = 0;
    graph_vertex **a_srknodes, **b_srknodes;
    graph_vertex **a_special = NULL, **b_special = NULL;
    graph_vertex *v;
    int acount, bcount, ascount, bscount;
    graph_ghtree_node *a_cut, *b_cut;
    graph_ghtree_node *parent, *child, *c, *nextc;
    double cutvalue;
    graph_ghtree_node *newnode;
    graph_vertex *anode     = NULL;
    graph_vertex *bnode     = NULL;
    int *elist              = NULL;
    double *ecap            = NULL;
    graph_vertex **cutverts = NULL;
    int i, ecount, cutcount;
    double amaxval, bmaxval;
    graph_vertex *amaxv = NULL, *bmaxv = NULL, *oldmax;
    int oldmcount;

    /***********************************************************/
    /* termination cases for the recursion                     */
    /***********************************************************/
    if (scount < 2)
        return 0;

    /***********************************************************/
    /*Select s and t randomly                                  */
    /***********************************************************/
    anode = special[rand() % scount];
    bnode = special[rand() % scount];
    while (anode->i == bnode->i) bnode = special[rand() % scount];

    /***********************************************************/
    /* Solve (s-t)-minimum cut*/
    /***********************************************************/
    rval =
    graph_get_mincut_st(graph, anode, bnode, &cutvalue, &cutverts, &cutcount);
    check_rval(rval, " failed", CLEANUP);

    /***********************************************************/
    /* Divide the search in two                                */
    /***********************************************************/

    /* mark the cut */
    graph->marker++;
    int j;
    for (i = 0; i < cutcount; i++)
    {
        cutverts[i]->mark = graph->marker;
        for (j = 1; j < ghtree->n[cutverts[i]->i]->mcount; j++)
            ghtree->n[cutverts[i]->i]->members[j]->special->mark =
            graph->marker;
    }
    for (i = 0; i < scount; i++) special[i]->mark_aux = graph->marker;

    /* divide them up */
    a_srknodes = malloc((graph->n3v - cutcount) * sizeof(graph_vertex *));
    a_special  = malloc((graph->n3v - cutcount) * sizeof(graph_vertex *));
    b_srknodes = malloc(cutcount * sizeof(graph_vertex *));
    b_special  = malloc(cutcount * sizeof(graph_vertex *));
    acount = bcount = 0;
    ascount = bscount = 0;
    amaxval           = -SOLVER_MAXDOUBLE;
    bmaxval           = -SOLVER_MAXDOUBLE;
    for (v = graph->tail; v; v = v->next)
    {
        if (v->mark_aux == graph->marker && v->mark != graph->marker)
        {
            a_special[ascount++] = v;
            a_srknodes[acount++] = v;
            if (amaxval < v->max->orig->y)
            {
                amaxval = v->max->orig->y;
                amaxv   = v->max;
            }
        }
        else if (v->mark_aux == graph->marker && v->mark == graph->marker)
        {
            b_special[bscount++] = v;
            b_srknodes[bcount++] = v;
            if (bmaxval < v->max->orig->y)
            {
                bmaxval = v->max->orig->y;
                bmaxv   = v->max;
            }
        }
        else if (v->mark_aux < graph->marker && v->mark != graph->marker)
        {
            a_srknodes[acount++] = v;
            if (amaxval < v->max->orig->y)
            {
                amaxval = v->max->orig->y;
                amaxv   = v->max;
            }
        }
        else if (v->mark_aux < graph->marker && v->mark == graph->marker)
        {
            b_srknodes[bcount++] = v;
            if (bmaxval < v->max->orig->y)
            {
                bmaxval = v->max->orig->y;
                bmaxv   = v->max;
            }
        }
        else
        {
            exit(1);
        }
    }

    assert(acount + bcount == graph->n3v);

    if (!ascount)
    {
        fprintf(stderr, "a_special is null\n");
        if (!bscount)
            fprintf(stderr, "And so is b_special\n");
        rval = 1;
        goto CLEANUP;
    }
    else if (!bscount)
    {
        fprintf(stderr, "b_special is null\n");
        rval = 1;
        goto CLEANUP;
    }

    /***********************************************************/
    /* Reorder the ghtree nodes (maintain the root at the top) */
    /***********************************************************/
    if (n->special->mark != graph->marker)
    {
        newnode = ghtree->n[bnode->i];
        a_cut   = n;
        b_cut   = ghtree->n[bnode->i];
        if (!n->parent || n->parent->special->mark != graph->marker)
        {
            parent         = a_cut;
            child          = b_cut;
            parent->in_max = parent->special->max;
            child->in_max  = bmaxv;
        }
        else
        {
            parent = b_cut;
            child  = a_cut;

            if (child->next_sibling)
                child->next_sibling->prev_sibling = child->prev_sibling;
            if (child->prev_sibling)
                child->prev_sibling->next_sibling = child->next_sibling;
            else
                child->parent->child = child->next_sibling;

            parent->next_sibling = child->parent->child;
            if (child->parent->child)
                child->parent->child->prev_sibling = parent;
            child->parent->child = parent;
            parent->prev_sibling = NULL;
            parent->parent       = child->parent;

            child->next_sibling = NULL;
            child->prev_sibling = NULL;
            child->parent       = NULL;

            parent->cutval = child->cutval;
            child->cutval  = 0.0;
            child->in_max  = amaxv;
            parent->in_max = parent->special->max;
        }
    }
    else
    {
        newnode = ghtree->n[anode->i];
        a_cut   = ghtree->n[anode->i];
        b_cut   = n;
        if (!n->parent || n->parent->special->mark == graph->marker)
        {
            parent         = b_cut;
            child          = a_cut;
            child->in_max  = amaxv;
            parent->in_max = parent->special->max;
        }
        else
        {
            parent = a_cut;
            child  = b_cut;

            if (child->next_sibling)
                child->next_sibling->prev_sibling = child->prev_sibling;
            if (child->prev_sibling)
                child->prev_sibling->next_sibling = child->next_sibling;
            else
                child->parent->child = child->next_sibling;

            parent->next_sibling = child->parent->child;
            if (child->parent->child)
                child->parent->child->prev_sibling = parent;
            child->parent->child = parent;
            parent->prev_sibling = NULL;
            parent->parent       = child->parent;

            child->next_sibling = NULL;
            child->prev_sibling = NULL;
            child->parent       = NULL;

            parent->cutval = child->cutval;
            child->cutval  = 0.0;
            parent->in_max = parent->special->max;
            child->in_max  = bmaxv;
        }
    }

    /* divide the children up between the two sides */
    for (c = n->child, n->child = NULL, n->nchild = 0; c; c = nextc)
    {
        nextc = c->next_sibling;
        if ((c->special->mark == graph->marker &&
             n->special->mark == graph->marker) ||
            (c->special->mark != graph->marker &&
             n->special->mark != graph->marker))
        {
            c->next_sibling = n->child;
            if (n->child)
                n->child->prev_sibling = c;
            c->prev_sibling = NULL;
            n->child        = c;
            c->parent       = n;
            n->nchild++;
            if (n->in_max->orig->y < c->in_max->orig->y)
                n->in_max = c->in_max;
        }
        else
        {
            c->next_sibling = newnode->child;
            if (newnode->child)
                newnode->child->prev_sibling = c;
            c->prev_sibling = NULL;
            newnode->child  = c;
            newnode->nchild++;
            c->parent = newnode;
            if (newnode->in_max->orig->y < c->in_max->orig->y)
                newnode->in_max = c->in_max;
        }
    }

    child->parent = parent;
    child->cutval = cutvalue;
    if (parent->child)
        parent->child->prev_sibling = child;
    child->next_sibling = parent->child;
    child->prev_sibling = NULL;
    parent->child       = child;
    parent->nchild++;
    if (parent->in_max->orig->y < child->in_max->orig->y)
        parent->in_max = child->in_max;

    elist = malloc(2 * graph->na * sizeof(int));
    ecap  = malloc(graph->na * sizeof(double));

    /***********************************************************/
    /* A                                                       */
    /***********************************************************/
    oldmax    = b_cut->special->max;
    oldmcount = b_cut->mcount;
    rval      = shrinkdown(ghtree, graph, bcount, b_srknodes, b_cut->special,
                           &ecount, elist, ecap);
    check_rval(rval, "shrinkdown failed", CLEANUP);
    if (b_srknodes)
        free(b_srknodes);

    rval = gh_recursive(graph, ghtree, a_cut, ascount, a_special);
    if (rval)
        goto CLEANUP;

    unshrink(graph, b_cut->special, ecount, elist, ecap);
    b_cut->special->max = oldmax;
    b_cut->mcount       = oldmcount;

    /***********************************************************/
    /* B                                                       */
    /***********************************************************/
    oldmax    = a_cut->special->max;
    oldmcount = a_cut->mcount;
    rval      = shrinkdown(ghtree, graph, acount, a_srknodes, a_cut->special,
                           &ecount, elist, ecap);
    check_rval(rval, "shrinkdown failed", CLEANUP);
    if (a_srknodes)
        free(a_srknodes);

    rval = gh_recursive(graph, ghtree, b_cut, bscount, b_special);
    if (rval)
        goto CLEANUP;

    unshrink(graph, a_cut->special, ecount, elist, ecap);
    a_cut->special->max = oldmax;
    a_cut->mcount       = oldmcount;

    rval = 0;

CLEANUP:

    if (elist)
        free(elist);
    if (a_special)
        free(a_special);
    if (b_special)
        free(b_special);
    if (ecap)
        free(ecap);
    if (cutverts)
        free(cutverts);

    return rval;
}

graph_ghtree *
graph_create_ghtree(void)
{
    graph_ghtree *ghtree = malloc(sizeof(graph_ghtree));
    if (!ghtree)
    {
        printf("out of memory in create_ghtree");
        return NULL;
    }
    ghtree->nn   = 0;
    ghtree->n    = NULL;
    ghtree->root = NULL;
    return ghtree;
}

graph_ghtree_node *
graph_add_ghtree_node(graph_ghtree *ghtree)
{
    graph_ghtree_node *node;
    assert(0 <= ghtree->nn && ghtree->nn <= ghtree->nn);

    node = ghtree->n[ghtree->nn] = malloc(sizeof(graph_ghtree_node));
    node->i                      = ghtree->nn++;
    node->parent                 = NULL;
    node->next_sibling           = NULL;
    node->child                  = NULL;
    node->cutval                 = 0.0;
    return node;
}

void
graph_add_ghtree_nodes(graph_ghtree *ghtree, int nn)
{
    int i;
    graph_ghtree_node *node;
    ghtree->n = malloc(nn * sizeof(graph_ghtree_node *));

    for (i = 0; i < nn; i++)
    {
        node = ghtree->n[i] = malloc(sizeof(graph_ghtree_node));
        node->i             = ghtree->nn++;
        node->in_max        = NULL;
        node->parent        = NULL;
        node->next_sibling  = NULL;
        node->prev_sibling  = NULL;
        node->child         = NULL;
        node->cutval        = 0.0;
        node->nchild        = 0;
        node->ndesc         = 0;
        node->members       = malloc(nn * sizeof(graph_ghtree_node *));
        node->members[0]    = node;
        node->mcount        = 1;
    }
    ghtree->root = ghtree->n[0];

    ghtree->nn = nn;
}

void
graph_free_ghtree(graph_ghtree *ghtree)
{
    int i;
    if (ghtree)
    {
        for (i = 0; i < ghtree->nn; i++)
        {
            ghtree_free_node(ghtree->n[i]);
        }
        if (ghtree->n)
            free(ghtree->n);
        free(ghtree);
    }
}

static void
ghtree_free_node(graph_ghtree_node *n)
{
    if (n)
    {
        free(n->members);
        free(n);
    }
}

static int
countdescendants(graph_ghtree_node *n)
{
    graph_ghtree_node *p;
    int i;

    i = 1;
    for (p = n->child; p; p = p->next_sibling) i += countdescendants(p);
    n->ndesc = i;
    return i;
}

/* shrink everything in the list to srknode. */
static int
shrinkdown(graph_ghtree *ghtree, solver_graph *graph, int nodecount,
           graph_vertex **nodeset, graph_vertex *srknode, int *ecount,
           int *elist, double *ecap)
{
    int rval = 0;
    int i;
    graph_arc *e, *next;
    graph_vertex *v, *other;
    double x;
    int j;

    *ecount = 0;

    graph->marker++;
    for (i = 0; i < nodecount; i++) nodeset[i]->mark = graph->marker;

    for (e = (srknode)->edge; e; e = next)
    {
        next  = outnext(e, (srknode));
        other = otherend(e, (srknode));

        elist[2 * (*ecount)]     = (srknode)->i;
        elist[2 * (*ecount) + 1] = other->i;
        ecap[(*ecount)++]        = e->x;

        if (other->mark == graph->marker)
            graph_del_arc(graph, &e);
    }
    for (i = 0; i < nodecount; i++)
    {
        v = nodeset[i];
        if (v->i != srknode->i)
        {
            x = 0.0;
            if (v->i != srknode->i)
            {
                if (v->max->orig->y > srknode->max->orig->y)
                    srknode->max = v->max;
                for (j = 0; j < ghtree->n[v->i]->mcount; j++)
                    ghtree->n[srknode->i]
                    ->members[ghtree->n[srknode->i]->mcount + j] =
                    ghtree->n[v->i]->members[j];
                ghtree->n[srknode->i]->mcount += ghtree->n[v->i]->mcount;
                for (e = v->edge; e; e = next)
                {
                    next  = outnext(e, v);
                    other = otherend(e, v);
                    if (other->mark != graph->marker)
                    {
                        x               = e->x;
                        ecap[(*ecount)] = e->x;

                        elist[2 * (*ecount)]       = v->i;
                        elist[2 * (*ecount)++ + 1] = other->i;
                        graph_del_arc(graph, &e);

                        e = graph_add_arc(graph, other->i, srknode->i);
                        check_null(e, "add_srk_arc failed", CLEANUP);
                        e->tail->y += x / 2.0;
                        e->head->y += x / 2.0;
                        e->x += x;
                    }
                    else
                    {
                        elist[2 * (*ecount)]     = v->i;
                        elist[2 * (*ecount) + 1] = other->i;
                        ecap[(*ecount)++]        = e->x;
                        graph_del_arc(graph, &e);
                    }
                }
            }
        }
    }
CLEANUP:
    return rval;
}

static void
unshrink(solver_graph *graph, graph_vertex *srknode, int ecount, int *elist,
         double *ecap)
{
    int i;
    graph_arc *e, *next;

    /* delete the new edges we added */
    for (e = srknode->edge; e; e = next)
    {
        next = outnext(e, srknode);
        graph_del_arc(graph, &e);
    }

    /* put back the edges we deleted when shrinking */
    for (i = 0; i < ecount; i++)
    {
        e    = graph_add_arc(graph, elist[2 * i], elist[2 * i + 1]);
        e->x = ecap[i];
        e->tail->y += ecap[i] / 2.0;
        e->head->y += ecap[i] / 2.0;
    }
}

static void
update_outmax(graph_ghtree_node *node)
{
    graph_ghtree_node *child;

    if (node->parent)
    {
        node->out_max = node->parent->special;
        if (node->out_max->y < node->parent->out_max->y)
            node->out_max = node->parent->out_max;
    }
    else
        node->out_max = node->special;

    if (node->parent)
    {
        for (child = node->parent->child; child; child = child->next_sibling)
        {
            if (child->i != node->i)
            {
                if (node->out_max->y < child->in_max->y)
                    node->out_max = child->in_max;
            }
        }
    }

    for (child = node->child; child; child = child->next_sibling)
        update_outmax(child);
}

void
graph_print_ghtree(graph_ghtree *ghtree)
{
    if (ghtree)
    {
        printf("\nGOMORY-HU TREE\n");
        ghtree_print_work(ghtree->root);
    }
    printf("\n");
}

static void
ghtree_print_work(graph_ghtree_node *root)
{
    graph_ghtree_node *n;

    printf("T-%d-%d (%.2f): ", root->i, root->special->i, root->special->y);
    if (root->parent)
        printf("Parent %d, ", root->parent->i);
    else
        printf("Parent NULL, ");
    if (root->child)
        printf("Child %d, ", root->child->i);
    else
        printf("Child NULL, ");
    if (root->next_sibling)
        printf("N Sibling %d, ", root->next_sibling->i);
    else
        printf("N Sibling NULL, ");
    if (root->prev_sibling)
        printf("P Sibling %d, ", root->prev_sibling->i);
    else
        printf("P Sibling NULL, ");
    if (root->parent)
        printf("Cnt %d, Nchild %d Val %.2f, ", root->ndesc, root->nchild,
               root->cutval);
    else
        printf("Cnt %d, Nchild %d ROOT, ", root->ndesc, root->nchild);
    printf("OutMax %d (%.2f), InMax %d (%.2f)", root->out_max->i,
           root->out_max->y, root->in_max->i, root->in_max->y);
    printf("\n");

    for (n = root->child; n; n = n->next_sibling) ghtree_print_work(n);
}
