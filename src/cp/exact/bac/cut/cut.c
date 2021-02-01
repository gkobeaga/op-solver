#include "cp/cp.h"
#include "cp/exact/bac/bac.h"

static void
create_cut_work(cp_cut *cut)
{
    cut->hcount       = 0;
    cut->tcount       = 0;
    cut->rhs          = 0.0;
    cut->sense        = 'X';
    cut->handles      = NULL;
    cut->teeth        = NULL;
    cut->cliques      = NULL;
    cut->verts        = NULL;
    cut->vycoef       = 0.0;
    cut->next         = NULL;
    cut->prev         = NULL;
    cut->age          = 0;
    cut->branch       = 0;
    cut->min_depth    = -1;
    cut->handle_cid   = NULL;
    cut->teeth_cid    = NULL;
    cut->logical      = NULL;
    cut->cover_edge   = NULL;
    cut->cover_vertex = NULL;
    cut->connect      = NULL;
    cut->path         = NULL;
    cut->hash_next    = NULL;
    cut->hash_prev    = NULL;
    cut->skel         = NULL;
}

cp_cut *
cp_create_cut(void)
{
    cp_cut *cut = malloc(sizeof(cp_cut));
    create_cut_work(cut);
    return cut;
}

void
cp_free_cut(cp_cut **cut)
{
    int i;

    if (*cut)
    {
        if ((*cut)->hcount + (*cut)->tcount)
        {

            if ((*cut)->handle_cid)
                free((*cut)->handle_cid);
            if ((*cut)->teeth_cid)
                free((*cut)->teeth_cid);
            if ((*cut)->handles)
            {
                for (i = 0; i < (*cut)->hcount; i++)
                {
                    if ((*cut)->handles[i])
                        clique_free(&(*cut)->handles[i]);
                    free((*cut)->handles[i]);
                }
                free((*cut)->handles);
            }
            if ((*cut)->teeth)
            {
                for (i = 0; i < (*cut)->tcount; i++)
                {
                    if ((*cut)->teeth[i])
                        clique_free(&(*cut)->teeth[i]);
                    free((*cut)->teeth[i]);
                }
                free((*cut)->teeth);
            }
            if ((*cut)->verts)
                free((*cut)->verts);
            if ((*cut)->skel)
            {
                cp_free_cut_skeleton(&(*cut)->skel);
                free((*cut)->skel);
            }
        }
        else if ((*cut)->logical)
        {
            free((*cut)->logical);
        }
        else if ((*cut)->cover_edge)
        {
            free((*cut)->cover_edge->arcs);
            free((*cut)->cover_edge);
        }
        else if ((*cut)->cover_vertex)
        {
            clique_free(&(*cut)->cover_vertex->verts);
            free((*cut)->cover_vertex);
        }
        else if ((*cut)->connect)
        {
            clique_free(&(*cut)->connect->verts);
            free((*cut)->connect->verts);
            free((*cut)->connect);
        }
        else if ((*cut)->path)
        {
            free((*cut)->path->arcs);
            if ((*cut)->path->farcs)
                free((*cut)->path->farcs);
            free((*cut)->path);
        }
        free(*cut);
        *cut = NULL;
    }
}

double
cp_eval_cut(solver_graph *graph, cp_cut *cut)
{
    int i;
    double slack = 0.0;
    graph_arc *arc;
    graph_vertex *v;

    if (cut->tcount + cut->hcount)
    {
        slack -= cut->rhs;

        int nzcnt          = 0;
        graph_arc **nzlist = malloc(graph->na * sizeof(graph_arc *));
        cp_get_cut_arcs(graph, cut, &nzcnt, nzlist);

        for (i = 0; i < nzcnt; i++)
        {
            arc = nzlist[i];
            if (arc->coef)
            {
                slack += arc->coef * arc->x;
                arc->coef = 0;
            }
        }

        for (i = 0; i < 2 * cut->tcount; i++)
        {
            if (cut->verts[i] >= 0)
                slack += cut->vycoef * graph->v[cut->verts[i]]->y;
        }

        if (cut->sense == 'L')
            slack = -slack;

        free(nzlist);
    }
    else if (cut->logical)
    {
        arc = graph_find_arc_hash(graph->archash, cut->logical->arc[0],
                                  cut->logical->arc[1]);
        if (arc)
        {
            slack = arc->x - graph->v[cut->logical->v]->y;
            slack = -slack;
        }
    }
    else if (cut->cover_edge)
    {
        slack = -cut->rhs;
        graph->marker++;
        for (i = 0; i < cut->cover_edge->na; i++)
        {
            arc =
            graph_find_arc_hash(graph->archash, cut->cover_edge->arcs[2 * i],
                                cut->cover_edge->arcs[2 * i + 1]);
            if (arc)
                slack += arc->x;

            if (cut->cover_edge->strong)
            {
                v = graph->v[cut->cover_edge->arcs[2 * i]];
                if (v->mark != graph->marker)
                {
                    slack -= v->y;
                    v->mark = graph->marker;
                }

                v = graph->v[cut->cover_edge->arcs[2 * i + 1]];
                if (v->mark != graph->marker)
                {
                    slack -= v->y;
                    v->mark = graph->marker;
                }
            }
        }
        slack = -slack;
    }
    else if (cut->cover_vertex)
    {
        int j;
        slack = -cut->rhs;
        FOREACH_NODE_IN_CLIQUE (i, cut->cover_vertex->verts, j)
        {
            v = graph->v[i];
            slack += v->y;
        }
        slack = -slack;
    }
    else if (cut->path)
    {
        slack = 0.0;
        for (i = 0; i < cut->path->na; i++)
        {
            arc = graph_find_arc_hash(graph->archash, cut->path->arcs[2 * i],
                                      cut->path->arcs[2 * i + 1]);
            if (arc)
                slack += arc->x;
        }
        for (i = 0; i < cut->path->fna; i++)
        {
            arc = graph_find_arc_hash(graph->archash, cut->path->farcs[2 * i],
                                      cut->path->farcs[2 * i + 1]);
            if (arc)
                slack -= arc->x;
        }

        graph->marker++;
        for (i = 0; i < 2 * cut->path->na; i++)
        {
            if (graph->v[cut->path->arcs[i]]->mark == graph->marker)
                slack -= graph->v[cut->path->arcs[i]]->y;
            graph->v[cut->path->arcs[i]]->mark = graph->marker;
        }
        slack = -slack;
    }

    return slack;
}

void
cp_get_cut_arcs(solver_graph *graph, cp_cut *cut, int *nzcnt,
                graph_arc **nzlist)
{
    int marker;
    int i, tmp, k;
    graph_arc *arc;
    graph_vertex *v, *other;
    graph_clique *tooth, *handle;

    *nzcnt = 0;

    if (cut->tcount + cut->hcount)
    {
        for (i = 0; i < cut->hcount; i++)
        {
            if (cut->handles)
                handle = cut->handles[i];
            else
                handle = cut->cliques[cut->handle_cid[i]];
            (graph->marker)++;
            marker = graph->marker;

            FOREACH_NODE_IN_CLIQUE (k, handle, tmp)
            {
                graph->v[k]->mark = marker;
            }

            FOREACH_NODE_IN_CLIQUE (k, handle, tmp)
            {
                v = graph->v[k];
                if (graph->directed)
                {
                    for (arc = v->out; arc; arc = arc->t_next)
                    {
                        if (arc->head->mark != marker)
                        {
                            if (!arc->coef)
                                nzlist[(*nzcnt)++] = arc;
                            arc->coef += 1;
                        }
                    }
                    for (arc = v->in; arc; arc = arc->h_next)
                    {
                        if (arc->tail->mark != marker)
                        {
                            if (!arc->coef)
                                nzlist[(*nzcnt)++] = arc;
                            arc->coef += 1;
                        }
                    }
                }
                else
                {
                    for (arc = v->edge; arc; arc = outnext(arc, v))
                    {
                        other = otherend(arc, v);
                        if (other->mark != marker)
                        {
                            if (!arc->coef)
                                nzlist[(*nzcnt)++] = arc;
                            arc->coef += 1;
                        }
                    }
                }
            }
        }

        for (i = 0; i < cut->tcount; i++)
        {
            if (cut->teeth)
                tooth = cut->teeth[i];
            else
                tooth = cut->cliques[cut->teeth_cid[i]];
            (graph->marker)++;
            marker = graph->marker;
            FOREACH_NODE_IN_CLIQUE (k, tooth, tmp)
            {
                graph->v[k]->mark = marker;
            }

            FOREACH_NODE_IN_CLIQUE (k, tooth, tmp)
            {
                v = graph->v[k];
                if (graph->directed)
                {
                    for (arc = v->out; arc; arc = arc->t_next)
                    {
                        if (arc->head->mark != marker)
                        {
                            if (!arc->coef)
                                nzlist[(*nzcnt)++] = arc;
                            arc->coef += 1;
                        }
                    }
                    for (arc = v->in; arc; arc = arc->h_next)
                    {
                        if (arc->tail->mark != marker)
                        {
                            if (!arc->coef)
                                nzlist[(*nzcnt)++] = arc;
                            arc->coef += 1;
                        }
                    }
                }
                else
                {
                    for (arc = v->edge; arc; arc = outnext(arc, v))
                    {
                        other = otherend(arc, v);
                        if (other->mark != marker)
                        {
                            if (!arc->coef)
                                nzlist[(*nzcnt)++] = arc;
                            arc->coef += 1;
                        }
                    }
                }
            }
        }
    }
    else if (cut->logical)
    {
        arc = graph_find_arc_hash(graph->archash, cut->logical->arc[0],
                                  cut->logical->arc[1]);
        if (arc)
        {
            if (!arc->coef)
                nzlist[(*nzcnt)++] = arc;
            arc->coef += 1;
        }
    }
    else if (cut->cover_edge)
    {
        for (i = 0; i < cut->cover_edge->na; i++)
        {
            arc =
            graph_find_arc_hash(graph->archash, cut->cover_edge->arcs[2 * i],
                                cut->cover_edge->arcs[2 * i + 1]);
            if (arc)
            {
                if (!arc->coef)
                    nzlist[(*nzcnt)++] = arc;
                arc->coef += 1;
            }
        }
    }
    else if (cut->path)
    {
        for (i = 0; i < cut->path->na; i++)
        {
            arc = graph_find_arc_hash(graph->archash, cut->path->arcs[2 * i],
                                      cut->path->arcs[2 * i + 1]);
            if (arc)
            {
                if (!arc->coef)
                    nzlist[(*nzcnt)++] = arc;
                arc->coef += 1;
            }
        }
        for (i = 0; i < cut->path->fna; i++)
        {
            arc = graph_find_arc_hash(graph->archash, cut->path->farcs[2 * i],
                                      cut->path->farcs[2 * i + 1]);
            if (arc)
            {
                if (!arc->coef)
                    nzlist[(*nzcnt)++] = arc;
                arc->coef -= 1;
            }
        }
    }
}

int
cp_add_lp_cut(cp_prob *prob, cp_exact_bac_env *bac_env, cp_cut *cut)
{
    cp_cut_repo *cuts = bac_env->cuts;

    if (cut->hcount + cut->tcount)
    {
        cp_register_cut_repo_cliques(cuts, cut);
        if (cp_find_cut_hash(cuts->hash, cut))
        {
            cp_unregister_cut_repo_cliques(cuts, cut);
            if (cut->skel)
            {
                cp_free_cut_skeleton(&cut->skel);
            }
            return 0;
        }
        cp_add_cut_hash(cuts->hash, cut);
    }

    if (cuts->count + 1 >= cuts->space)
    {
        realloc_scale(cuts->cuts, cuts->space, cuts->count + 1, 1.3);
    }
    cuts->cuts[cuts->count++] = cut;

    return cuts->count;
}

int
cp_del_lp_cut(cp_prob *prob, cp_exact_bac_env *bac_env, int ind)
{
    int i;

    cp_unregister_cut_repo_cliques(bac_env->cuts, bac_env->cuts->cuts[ind]);
    cp_free_cut(&(bac_env->cuts->cuts[ind]));

    for (i = ind + 1; i < bac_env->cuts->count; i++)
        bac_env->cuts->cuts[i - 1] = bac_env->cuts->cuts[i];

    bac_env->cuts->count--;
    return 0;
}

void
cp_print_cut(cp_cut *cut)
{
    int i;
    if (cut->hcount + cut->tcount)
    {
        printf("CP Cut:\n");
        if (cut->tcount)
        {
            printf(" - Verts:");
            for (i = 0; i < cut->tcount; i++)
            {
                printf("  %d, %d\n", cut->verts[2 * i], cut->verts[2 * i + 1]);
            }
        }
        if (cut->hcount)
        {
            printf(" - Handle: ");
            if (cut->handles)
                clique_print(cut->handles[0]);
            else
                clique_print(cut->cliques[cut->handle_cid[0]]);
        }
        if (cut->tcount)
        {
            printf(" - Teeth [%d]\n", cut->tcount);
            for (i = 0; i < cut->tcount; i++)
            {
                printf("   + Tooth-%d: ", i);
                if (cut->teeth)
                    clique_print(cut->teeth[i]);
                else
                    clique_print(cut->cliques[cut->teeth_cid[i]]);
            }
        }
    }
    else if (cut->logical)
    {
        printf("\nCP Cut:\n");
        printf(" - Logical: ");
        printf(" - Arc: (%d, %d)\n", cut->logical->arc[0],
               cut->logical->arc[1]);
        printf(" - Vert: %d\n\n", cut->logical->v);
    }
}

cp_cut *
cp_conv_cut_sol2connect(cp_sol *sol)
{
    int rval = 0;

    if (sol->ns == 0)
        return NULL;

    cp_cut *cut = cp_create_cut();
    check_null(cut, "out of memory", cleanup);

    cut->hcount = 1;
    cut->tcount = 0;

    cut->handles = malloc(sizeof(graph_clique *));
    check_null(cut->handles, "out of memory", cleanup);
    rval = clique_conv_array2clique(sol->cycle, sol->ns, &cut->handles[0]);
    check_rval(rval, "failed", cleanup);

    cut->rhs    = 2.0;
    cut->sense  = 'G';
    cut->branch = 0;

    rval = cp_build_cut_skeleton(cut, sol->tot_n);
    check_rval(rval, "failed", cleanup);

    return cut;

cleanup:

    if (rval)
    {
        if (cut)
        {
            cp_free_cut(&cut);
        }
    }

    return NULL;
}

cp_cut *
cp_conv_cut_verts2connect(solver_graph *graph, int vcount, graph_vertex **verts)
{
    int rval    = 0;
    cp_cut *cut = cp_create_cut();
    check_null(cut, "out of memory", cleanup);

    cut->hcount = 1;
    cut->tcount = 0;

    cut->handles = malloc(sizeof(graph_clique *));
    check_null(cut->handles, "out of memory", cleanup);
    cut->handles[0] = clique_conv_vertices2clique(graph, verts, vcount);

    assert(cut->handles[0]->nodes[0].lo == 0);

    cut->rhs    = 2.0;
    cut->sense  = 'G';
    cut->branch = 0;

    rval = cp_build_cut_skeleton(cut, graph->nv);
    check_rval(rval, "", cleanup);

    return cut;

cleanup:

    if (rval && cut)
        cp_free_cut(&cut);

    return NULL;
}

cp_cut_list *
cp_create_cut_list(int count, double max_val)
{
    int rval = 0;
    int i;
    cp_cut_list *list = NULL;

    list = malloc(sizeof(cp_cut_list));
    check_null(list, "out of memory", cleanup);

    list->cuts = malloc((count + 1) * sizeof(cp_cut_list_item));
    check_null(list, "out of memory", cleanup);

    for (i = 0; i < count; i++)
    {
        list->cuts[i].val = 0.0;
        list->cuts[i].cut = NULL;
    }

    list->cuts[count].val = SOLVER_MAXDOUBLE;
    list->cuts[count].cut = NULL;
    list->max_val         = max_val;
    list->max_size        = count;

    return list;

cleanup:
    if (rval)
        fprintf(stderr, "create_cutlist failed\n");
    return NULL;
}

void
cp_insert_cut_list(cp_cut_list *list, double val, cp_cut *cut, int *finish)
{
    int k;

    *finish = 0;
    if (list->cuts[0].val < val)
    {
        cp_free_cut(&(list->cuts[0].cut));
        for (k = 0; list->cuts[k + 1].val < val; k++)
        {
            list->cuts[k].cut = list->cuts[k + 1].cut;
            list->cuts[k].val = list->cuts[k + 1].val;
        }
        list->cuts[k].val = val;
        list->cuts[k].cut = cut;
        if (list->cuts[0].val > list->max_val - SOLVER_ZEROPLUS)
            *finish = 1;
    }
    else
        cp_free_cut(&cut);
}

void
cp_get_cut_list(cp_prob *cp, cp_cut_list *cut_list, int *cutcount,
                cp_cut **cuts)
{
    int rval = 0;
    int i;
    for (i = cut_list->max_size - 1; i >= 0; i--)
    {
        if (cut_list->cuts[i].cut)
        {
            rval = cp_build_cut_skeleton(cut_list->cuts[i].cut, cp->n);
            check_rval(rval, "failed", cleanup);

            cut_list->cuts[i].cut->next = *cuts;
            *cuts                       = cut_list->cuts[i].cut;
            cut_list->cuts[i].cut       = NULL;
            (*cutcount)++;
        }
    }

    return;

cleanup:
    return;
}

void
cp_append_cut_list(cp_cut_list *main_list, cp_cut_list *from_list, int *stop)
{

    int ind = main_list->max_size - 1;
    while (ind >= 0 && main_list->cuts[ind].cut) ind--;
    if (ind >= 0)
    {
        for (int i = from_list->max_size - 1; i >= 0 && ind >= 0; i--)
        {
            if (from_list->cuts[i].cut)
            {
                main_list->cuts[ind].cut = from_list->cuts[i].cut;
                main_list->cuts[ind].val = from_list->cuts[i].val;
                ind--;
                from_list->cuts[i].cut = NULL;
                from_list->cuts[i].val = 0;
            }
        }
    }
    if (ind <= 0)
        *stop = 1;

    return;
}

void
cp_erase_cut_list(cp_cut_list *cut_list)
{
    for (int i = cut_list->max_size - 1; i >= 0; i--)
    {
        cut_list->cuts[i].cut = NULL;
        cut_list->cuts[i].val = 0.0;
    }

    return;
}

void
cp_free_cut_list(cp_cut_list **cut_list)
{
    if (*cut_list)
    {
        if ((*cut_list)->cuts)
            free((*cut_list)->cuts);
        free(*cut_list);
        *cut_list = NULL;
    }
}
