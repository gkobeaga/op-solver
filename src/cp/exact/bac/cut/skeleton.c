#include "cp/cp.h"
#include "cp/exact/bac/bac.h"

#define SKELETON_WILD 0

#undef DEBUG_CONSTRUCT

cp_cut_skeleton *
cp_create_cut_skeleton(void)
{
    cp_cut_skeleton *skel = malloc(sizeof(cp_cut_skeleton));
    skel->atomcount       = 0;
    skel->atoms           = NULL;
    skel->nverts          = 0;
    skel->verts           = NULL;
    return skel;
}

void
cp_free_cut_skeleton(cp_cut_skeleton **skel)
{
    if (*skel)
    {
        (*skel)->atomcount = 0;
        if ((*skel)->atoms)
            free((*skel)->atoms);
        (*skel)->nverts = 0;
        if ((*skel)->verts)
            free((*skel)->verts);
        free(*skel);
        *skel = NULL;
    }
}

int
cp_copy_cut_skeleton(cp_cut_skeleton *in, cp_cut_skeleton *out)
{
    int i;

    if (in->atomcount == 0)
        return 0;
    if (out->atoms)
        free(out->atoms);
    out->atoms = malloc(in->atomcount * sizeof(int));
    if (!out->atoms)
    {
        fprintf(stderr, "Out of memory\n");
        return 1;
    }

    for (i = 0; i < in->atomcount; i++) out->atoms[i] = in->atoms[i];

    out->atomcount = in->atomcount;
    if (in->nverts)
    {
        out->nverts = in->nverts;
        if (out->verts)
            free(out->verts);
        out->verts = malloc(out->nverts * sizeof(int));
        for (i = 0; i < in->nverts; i++) out->verts[i] = in->verts[i];
    }

    return 0;
}

static int
sort_int(const void *a, const void *b)
{
    return (*(int *)a - *(int *)b);
}

int
cp_build_cut_skeleton(cp_cut *cut, int nodecount)
{
    int rval   = 0;
    int *label = NULL;
    int ccount;
    int *cnodes = NULL;
    int atomcount;
    int atomcount_save;
    int *atomsize = NULL;
    int *atomnew  = NULL;
    int *atomwork = NULL;
    int *atoms    = NULL;
    int i;
    int j;
    int tmp;

    if (cut->hcount + cut->tcount == 0)
        return 0;

    if (!(cut->skel))
        cut->skel = cp_create_cut_skeleton();

    label = malloc(nodecount * sizeof(int));
    if (!label)
    {
        fprintf(stderr, "Out of memory\n");
        rval = 1;
        goto CLEANUP;
    }

    for (i = 0; i < cut->hcount; i++)
    {
        graph_clique *handle = cut->handles[i];
        FOREACH_NODE_IN_CLIQUE (j, handle, tmp)
        {
            label[j] = 0;
        }
    }

    for (i = 0; i < cut->tcount; i++)
    {
        graph_clique *tooth = cut->teeth[i];
        FOREACH_NODE_IN_CLIQUE (j, tooth, tmp)
        {
            label[j] = 0;
        }
    }

    ccount = 0;
    for (i = 0; i < cut->hcount; i++)
    {
        graph_clique *handle = cut->handles[i];
        FOREACH_NODE_IN_CLIQUE (j, handle, tmp)
        {
            if (label[j] == 0)
            {
                label[j] = 1;
                ccount++;
            }
        }
    }
    for (i = 0; i < cut->tcount; i++)
    {
        graph_clique *tooth = cut->teeth[i];
        FOREACH_NODE_IN_CLIQUE (j, tooth, tmp)
        {
            if (label[j] == 0)
            {
                label[j] = 1;
                ccount++;
            }
        }
    }

    cnodes   = malloc(ccount * sizeof(int));
    atomsize = malloc((ccount + 1) * sizeof(int));
    atomnew  = malloc((ccount + 1) * sizeof(int));
    atomwork = malloc((ccount + 1) * sizeof(int));
    if (!cnodes || !atomsize || !atomnew || !atomwork)
    {
        fprintf(stderr, "Out of memory\n");
        rval = 1;
        goto CLEANUP;
    }

    /* collect nodes */
    ccount = 0;
    for (i = 0; i < cut->hcount; i++)
    {
        graph_clique *handle = cut->handles[i];
        FOREACH_NODE_IN_CLIQUE (j, handle, tmp)
        {
            if (label[j] == 1)
            {
                label[j]         = 0;
                cnodes[ccount++] = j;
            }
        }
    }
    for (i = 0; i < cut->tcount; i++)
    {
        graph_clique *tooth = cut->teeth[i];
        FOREACH_NODE_IN_CLIQUE (j, tooth, tmp)
        {
            if (label[j] == 1)
            {
                label[j]         = 0;
                cnodes[ccount++] = j;
            }
        }
    }

    qsort(cnodes, ccount, sizeof(int), sort_int);

    /* refine atoms */
    atomsize[0] = ccount;
    atomcount   = 1;
    for (i = 0; i < cut->hcount; i++)
    {
        graph_clique *handle = cut->handles[i];
        for (j = 0; j < atomcount; j++) atomwork[j] = 0;

        FOREACH_NODE_IN_CLIQUE (j, handle, tmp)
        {
            atomwork[label[j]]++;
        }

        atomcount_save = atomcount;
        for (j = 0; j < atomcount_save; j++)
        {
            if (atomwork[j] == 0)
            {
                atomnew[j] = -1;
            }
            else if (atomwork[j] == atomsize[j])
            {
                atomnew[j] = j;
            }
            else
            {
                atomsize[atomcount] = atomwork[j];
                atomsize[j] -= atomwork[j];
                atomnew[j] = atomcount;
                atomcount++;
            }
        }
        FOREACH_NODE_IN_CLIQUE (j, handle, tmp)
        {
            label[j] = atomnew[label[j]];
        }
    }
    for (i = 0; i < cut->tcount; i++)
    {
        graph_clique *tooth = cut->teeth[i];
        for (j = 0; j < atomcount; j++) atomwork[j] = 0;

        FOREACH_NODE_IN_CLIQUE (j, tooth, tmp)
        {
            atomwork[label[j]]++;
        }

        atomcount_save = atomcount;
        for (j = 0; j < atomcount_save; j++)
        {
            if (atomwork[j] == 0)
            {
                atomnew[j] = -1;
            }
            else if (atomwork[j] == atomsize[j])
            {
                atomnew[j] = j;
            }
            else
            {
                atomsize[atomcount] = atomwork[j];
                atomsize[j] -= atomwork[j];
                atomnew[j] = atomcount;
                atomcount++;
            }
        }
        FOREACH_NODE_IN_CLIQUE (j, tooth, tmp)
        {
            label[j] = atomnew[label[j]];
        }
    }

    atomcount_save = atomcount;
    if (ccount < nodecount)
    {
        atomcount += 1;
    }

    atoms = malloc(atomcount * sizeof(int));
    check_null(atoms, "Out of memory", CLEANUP);

    for (i = 0; i < atomcount; i++) atoms[i] = -1;

    for (i = 0; i < ccount; i++)
    {
        if (atoms[label[cnodes[i]]] == -1)
            atoms[label[cnodes[i]]] = cnodes[i];
    }

    if (ccount < nodecount)
    {
        if (cnodes[ccount - 1] == ccount - 1)
        {
            atoms[atomcount_save] = ccount;
        }
        else
        {
            for (i = 0; i < ccount; i++)
            {
                if (cnodes[i] != i)
                {
                    atoms[atomcount_save] = i;
                    break;
                }
            }
        }
    }

    qsort(atoms, atomcount, sizeof(int), sort_int);

    cut->skel->atoms     = atoms;
    cut->skel->atomcount = atomcount;
    atoms                = NULL;

    if (cut->verts)
    {
        cut->skel->nverts = 2 * cut->tcount;
        cut->skel->verts  = malloc(cut->skel->nverts * sizeof(int));
        for (i = 0; i < 2 * cut->tcount; i++)
            cut->skel->verts[i] = cut->verts[i];

        qsort(cut->skel->verts, cut->skel->nverts, sizeof(int), sort_int);
    }

    rval = 0;

CLEANUP:
    if (atoms)
        free(atoms);
    if (atomwork)
        free(atomwork);
    if (atomnew)
        free(atomnew);
    if (atomsize)
        free(atomsize);
    if (cnodes)
        free(cnodes);
    if (label)
        free(label);
    if (rval)
        cp_free_cut_skeleton(&(cut->skel));
    return rval;
}

int
cp_eq_cut_skeletons(cp_cut_skeleton *a, cp_cut_skeleton *b)
{
    int i;

    if (a->atomcount != b->atomcount)
        return 0;
    if (a->nverts != b->nverts)
        return 0;
    for (i = 0; i < a->nverts; i++)
    {
        if (a->verts[i] != b->verts[i])
            return 0;
    }

    for (i = 0; i < a->atomcount; i++)
    {
        if (a->atoms[i] != b->atoms[i])
            return 0;
    }

    return 1;
}
