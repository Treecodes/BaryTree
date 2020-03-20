#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../run_params/struct_run_params.h"
#include "../particles/struct_particles.h"

#include "struct_tree_linked_list_node.h"
#include "tree_linked_list.h"

#include "struct_tree.h"
#include "tree.h"


void Tree_Sources_Construct(struct Tree **tree_addr, struct Particles *sources, struct RunParams *run_params)
{
    struct TreeLinkedListNode *tree_linked_list = NULL;
    
    double xyzminmax[6];
    int numnodes = 0;
    int numleaves = 0;
    
    xyzminmax[0] = minval(sources->x, sources->num);
    xyzminmax[1] = maxval(sources->x, sources->num);
    xyzminmax[2] = minval(sources->y, sources->num);
    xyzminmax[3] = maxval(sources->y, sources->num);
    xyzminmax[4] = minval(sources->z, sources->num);
    xyzminmax[5] = maxval(sources->z, sources->num);
    
    TreeLinkedList_Sources_Construct(&tree_linked_list, sources, 1, sources->num,
            run_params->max_per_source_leaf, xyzminmax, 0, &numnodes, &numleaves);
    
    TreeLinkedList_SetIndex(tree_linked_list, 0);
    
    Tree_Alloc(tree_addr, numnodes);
    Tree_Fill(*tree_addr, tree_linked_list);
    (*tree_addr)->numleaves = numleaves;
    
    TreeLinkedList_Free(tree_linked_list);

    return;
}



void Tree_Targets_Construct(struct Tree **tree_addr, struct Particles *targets, struct RunParams *run_params)
{
    struct TreeLinkedListNode *tree_linked_list = NULL;
    
    double xyzminmax[6];
    int numnodes = 0;
    int numleaves = 0;
    
    xyzminmax[0] = minval(targets->x, targets->num);
    xyzminmax[1] = maxval(targets->x, targets->num);
    xyzminmax[2] = minval(targets->y, targets->num);
    xyzminmax[3] = maxval(targets->y, targets->num);
    xyzminmax[4] = minval(targets->z, targets->num);
    xyzminmax[5] = maxval(targets->z, targets->num);
    
    TreeLinkedList_Targets_Construct(&tree_linked_list, targets, 1, targets->num,
                    run_params->max_per_source_leaf, xyzminmax, 0, &numnodes, &numleaves);
    
    TreeLinkedList_SetIndex(tree_linked_list, 0);
    
    Tree_Alloc(tree_addr, numnodes);
    Tree_Fill(*tree_addr, tree_linked_list);
    (*tree_addr)->numleaves = numleaves;
    
    TreeLinkedList_Free(tree_linked_list);

    return;
}



void Tree_Alloc(struct Tree **tree_addr, int length)
{
    *tree_addr = malloc(sizeof(struct Tree));
    struct Tree *tree = *tree_addr;

    tree->numnodes = length;
    make_vector(tree->ibeg, length);
    make_vector(tree->iend, length);
    make_vector(tree->numpar, length);
    make_vector(tree->x_mid, length);
    make_vector(tree->y_mid, length);
    make_vector(tree->z_mid, length);
    make_vector(tree->x_min, length);
    make_vector(tree->y_min, length);
    make_vector(tree->z_min, length);
    make_vector(tree->x_max, length);
    make_vector(tree->y_max, length);
    make_vector(tree->z_max, length);
    make_vector(tree->level, length);
    make_vector(tree->cluster_ind, length);
    make_vector(tree->radius, length);
    make_vector(tree->num_children, length);
    make_vector(tree->children, 8*length);
    
    return;
}   /* END of function allocate_tree */



void Tree_Free(struct Tree *tree)
{
    if (tree != NULL) {
        free_vector(tree->ibeg);
        free_vector(tree->iend);
        free_vector(tree->numpar);
        free_vector(tree->x_mid);
        free_vector(tree->y_mid);
        free_vector(tree->z_mid);
        free_vector(tree->x_min);
        free_vector(tree->y_min);
        free_vector(tree->z_min);
        free_vector(tree->x_max);
        free_vector(tree->y_max);
        free_vector(tree->z_max);
        free_vector(tree->level);
        free_vector(tree->cluster_ind);
        free_vector(tree->radius);
        free_vector(tree->num_children);
        free_vector(tree->children);
        free(tree);
    }

    tree = NULL;

    return;
}   /* END of function allocate_tree */



void Tree_Realloc(struct Tree *tree, int newlength)
{
    tree->numnodes = newlength;
    realloc_vector(tree->ibeg, newlength);
    realloc_vector(tree->iend, newlength);
    realloc_vector(tree->numpar, newlength);
    realloc_vector(tree->x_mid, newlength);
    realloc_vector(tree->y_mid, newlength);
    realloc_vector(tree->z_mid, newlength);
    realloc_vector(tree->x_min, newlength);
    realloc_vector(tree->y_min, newlength);
    realloc_vector(tree->z_min, newlength);
    realloc_vector(tree->x_max, newlength);
    realloc_vector(tree->y_max, newlength);
    realloc_vector(tree->z_max, newlength);
    realloc_vector(tree->level, newlength);
    realloc_vector(tree->cluster_ind, newlength);
    realloc_vector(tree->radius, newlength);

    realloc_vector(tree->num_children, newlength);
    realloc_vector(tree->children, 8*newlength);

    return;
}   /* END of function allocate_tree */



void Tree_Fill(struct Tree *tree, struct TreeLinkedListNode *p)
{
    tree->x_mid[p->node_index] = p->x_mid;
    tree->y_mid[p->node_index] = p->y_mid;
    tree->z_mid[p->node_index] = p->z_mid;
    
    tree->x_min[p->node_index] = p->x_min;
    tree->y_min[p->node_index] = p->y_min;
    tree->z_min[p->node_index] = p->z_min;
    
    tree->x_max[p->node_index] = p->x_max;
    tree->y_max[p->node_index] = p->y_max;
    tree->z_max[p->node_index] = p->z_max;
    
    tree->ibeg[p->node_index] = p->ibeg;
    tree->iend[p->node_index] = p->iend;
    tree->numpar[p->node_index] = p->numpar;
    tree->level[p->node_index] = p->level;
    tree->radius[p->node_index] = p->radius;
    tree->cluster_ind[p->node_index] = p->node_index;

    tree->num_children[p->node_index] = p->num_children;

    for (int i = 0; i < p->num_children; i++) {
        tree->children[8*p->node_index+i] = (p->child[i])->node_index;
        Tree_Fill(tree, p->child[i]);
    }
    
    return;
} /* END of function Tree_CreateArray */
