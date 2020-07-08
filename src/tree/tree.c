#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <limits.h>
#include <mpi.h>

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
    
    int min_leaf_size = INT_MAX;
    int max_leaf_size = 0;
    
    xyzminmax[0] = minval(sources->x, sources->num);
    xyzminmax[1] = maxval(sources->x, sources->num);
    xyzminmax[2] = minval(sources->y, sources->num);
    xyzminmax[3] = maxval(sources->y, sources->num);
    xyzminmax[4] = minval(sources->z, sources->num);
    xyzminmax[5] = maxval(sources->z, sources->num);
    
    TreeLinkedList_Sources_Construct(&tree_linked_list, sources, 1, sources->num,
                    run_params->max_per_source_leaf, xyzminmax, &numnodes, &numleaves,
                    &min_leaf_size, &max_leaf_size);
    
    TreeLinkedList_SetIndex(tree_linked_list, 0);
    
    Tree_Alloc(tree_addr, numnodes);
    Tree_Fill(*tree_addr, tree_linked_list);
    (*tree_addr)->numleaves = numleaves;
    
    (*tree_addr)->min_leaf_size = min_leaf_size;
    (*tree_addr)->max_leaf_size = max_leaf_size;
    
    TreeLinkedList_Free(&tree_linked_list);

    return;
}



void Tree_Targets_Construct(struct Tree **tree_addr, struct Particles *targets, struct RunParams *run_params)
{
    struct TreeLinkedListNode *tree_linked_list = NULL;
    
    double xyzminmax[6];
    int xyzdim[3], xyzind[6];

    int numnodes = 0;
    int numleaves = 0;
    
    int min_leaf_size = INT_MAX;
    int max_leaf_size = 0;
    
    xyzminmax[0] = targets->xmin;
    xyzminmax[1] = targets->xmax;
    xyzminmax[2] = targets->ymin;
    xyzminmax[3] = targets->ymax;
    xyzminmax[4] = targets->zmin;
    xyzminmax[5] = targets->zmax;

    xyzdim[0] = targets->xdim;
    xyzdim[1] = targets->ydim;
    xyzdim[2] = targets->zdim;
    
    xyzind[0] = 0;
    xyzind[1] = targets->xdim-1;
    xyzind[2] = 0;
    xyzind[3] = targets->ydim-1;
    xyzind[4] = 0;
    xyzind[5] = targets->zdim-1;
    
    TreeLinkedList_Targets_Construct(&tree_linked_list,
                    run_params->max_per_target_leaf, xyzminmax, xyzdim, xyzind,
                    &numnodes, &numleaves, &min_leaf_size, &max_leaf_size);
    
    TreeLinkedList_SetIndex(tree_linked_list, 0);
    
    Tree_Alloc(tree_addr, numnodes);
    Tree_Fill(*tree_addr, tree_linked_list);
    (*tree_addr)->numleaves = numleaves;
    
    (*tree_addr)->min_leaf_size = min_leaf_size;
    (*tree_addr)->max_leaf_size = max_leaf_size;
    
    TreeLinkedList_Free(&tree_linked_list);

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
    make_vector(tree->cluster_ind, length);
    make_vector(tree->used, length);
    make_vector(tree->radius, length);
    make_vector(tree->num_children, length);
    make_vector(tree->children, 8*length);

    make_vector(tree->x_dim, length);
    make_vector(tree->y_dim, length);
    make_vector(tree->z_dim, length);

    make_vector(tree->x_low_ind, length);
    make_vector(tree->y_low_ind, length);
    make_vector(tree->z_low_ind, length);

    make_vector(tree->x_high_ind, length);
    make_vector(tree->y_high_ind, length);
    make_vector(tree->z_high_ind, length);
    
    return;
}   /* END of function allocate_tree */



void Tree_Free(struct Tree **tree_addr)
{
    struct Tree *tree = *tree_addr;
    
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
        free_vector(tree->cluster_ind);
        free_vector(tree->used);
        free_vector(tree->radius);
        free_vector(tree->num_children);
        free_vector(tree->children);

        free_vector(tree->x_dim);
        free_vector(tree->y_dim);
        free_vector(tree->z_dim);

        free_vector(tree->x_low_ind);
        free_vector(tree->y_low_ind);
        free_vector(tree->z_low_ind);

        free_vector(tree->x_high_ind);
        free_vector(tree->y_high_ind);
        free_vector(tree->z_high_ind);

        free(tree);
    }

    tree = NULL;

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
    tree->radius[p->node_index] = p->radius;
    tree->cluster_ind[p->node_index] = p->node_index;

    tree->num_children[p->node_index] = p->num_children;

    tree->x_dim[p->node_index] = p->x_dim;
    tree->y_dim[p->node_index] = p->y_dim;
    tree->z_dim[p->node_index] = p->z_dim;

    tree->x_low_ind[p->node_index] = p->x_low_ind;
    tree->y_low_ind[p->node_index] = p->y_low_ind;
    tree->z_low_ind[p->node_index] = p->z_low_ind;

    tree->x_high_ind[p->node_index] = p->x_high_ind;
    tree->y_high_ind[p->node_index] = p->y_high_ind;
    tree->z_high_ind[p->node_index] = p->z_high_ind;

    for (int i = 0; i < p->num_children; i++) {
        tree->children[8*p->node_index+i] = (p->child[i])->node_index;
        Tree_Fill(tree, p->child[i]);
    }
    
    return;
} /* END of function Tree_CreateArray */



void Tree_Print(struct Tree *tree)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Barrier(MPI_COMM_WORLD);
   
    int global_num_nodes,  max_num_nodes,  min_num_nodes;
    int global_num_leaves, max_num_leaves, min_num_leaves;
    int global_min_leaf_size, global_max_leaf_size;

    MPI_Reduce(&(tree->numnodes),   &global_num_nodes, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(tree->numnodes),      &max_num_nodes, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(tree->numnodes),      &min_num_nodes, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&(tree->numleaves), &global_num_leaves, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(tree->numleaves),    &max_num_leaves, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(tree->numleaves),    &min_num_leaves, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&(tree->max_leaf_size),    &global_max_leaf_size, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&(tree->min_leaf_size),    &global_min_leaf_size, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    
    if (rank == 0) {
        printf("[BaryTree]\n");
        printf("[BaryTree] Tree information: \n");
        printf("[BaryTree]\n");
        printf("[BaryTree]          Cumulative tree nodes across all ranks: %d\n", global_num_nodes);
        printf("[BaryTree]             Maximum tree nodes across all ranks: %d\n", max_num_nodes);
        printf("[BaryTree]             Minimum tree nodes across all ranks: %d\n", min_num_nodes);
        printf("[BaryTree]                                           Ratio: %f\n",
               (double)max_num_nodes / (double)min_num_nodes);
        printf("[BaryTree]\n");
        printf("[BaryTree]         Cumulative tree leaves across all ranks: %d\n", global_num_leaves);
        printf("[BaryTree]            Maximum tree leaves across all ranks: %d\n", max_num_leaves);
        printf("[BaryTree]            Minimum tree leaves across all ranks: %d\n", min_num_leaves);
        printf("[BaryTree]                                           Ratio: %f\n",
               (double)max_num_leaves / (double)min_num_leaves);
        printf("[BaryTree]\n");
        printf("[BaryTree]         Maximum tree leaf size across all ranks: %d\n", global_max_leaf_size);
        printf("[BaryTree]         Minimum tree leaf size across all ranks: %d\n", global_min_leaf_size);
        printf("[BaryTree]                                           Ratio: %f\n",
               (double)global_max_leaf_size / (double)global_min_leaf_size);
        printf("[BaryTree]\n");
    }

    return;
}
