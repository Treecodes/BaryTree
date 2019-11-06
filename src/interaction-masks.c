#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "array.h"
#include "treedriver.h"
#include "tools.h"
#include "particles.h"
#include "sort.h"
#include "tnode.h"
#include "batch.h"
#include "tree.h"


void remote_interaction_lists(const struct tnode_array *tree_array, struct batch *batches,
					int *approx_list_unpacked,int *approx_list_packed, int *direct_list,
                    int *num_batch_approx, int *num_batch_direct)
{
    /* local variables */
	int rank, numProcs;

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);


    int i, j;

    int **batches_ind;
    double **batches_center;
    double *batches_radius;

    int tree_numnodes;
    const int *tree_numpar, *tree_level;
    const double *tree_x_mid, *tree_y_mid, *tree_z_mid, *tree_radius;
    const int *tree_num_children, *tree_children;

    batches_ind = batches->index;
    batches_center = batches->center;
    batches_radius = batches->radius;

    tree_numnodes = tree_array->numnodes;
    tree_numpar = tree_array->numpar;
    tree_level = tree_array->level;
    tree_radius = tree_array->radius;
    tree_x_mid = tree_array->x_mid;
    tree_y_mid = tree_array->y_mid;
    tree_z_mid = tree_array->z_mid;

    tree_num_children = tree_array->num_children;
    tree_children = tree_array->children;

    for (int i = 0; i < tree_numnodes; i++) approx_list_unpacked[i] = -1;
    for (int i = 0; i < tree_numnodes; i++) approx_list_packed[i] = -1;
   	for (int i = 0; i < tree_numnodes; i++) direct_list[i] = -1;


    int **temp_tree_inter_list, **temp_direct_inter_list;
    int *num_tree_inter, *num_direct_inter;
    int *sizeof_tree_inter_list, *sizeof_direct_inter_list;

    make_matrix(temp_tree_inter_list, batches->num, 50);
    make_matrix(temp_direct_inter_list, batches->num, 50);

    make_vector(num_tree_inter, batches->num);
    make_vector(num_direct_inter, batches->num);

    make_vector(sizeof_tree_inter_list, batches->num);
    make_vector(sizeof_direct_inter_list, batches->num);
   
    int loopsize = batches->num;
	for (int i = 0; i < loopsize; i++) num_tree_inter[i] = 0;
	for (int i = 0; i < loopsize; i++) num_direct_inter[i] = 0;

	for (int i = 0; i < loopsize; i++) sizeof_tree_inter_list[i] = 50;
	for (int i = 0; i < loopsize; i++) sizeof_direct_inter_list[i] = 50;
    
    for (int i = 0; i < loopsize; i++)
        for(int j = 0; j < 50; j++)
            temp_tree_inter_list[i][j] = -1;

    for (int i = 0; i < loopsize; i++)
        for(int j = 0; j < 50; j++)
            temp_direct_inter_list[i][j] = -1;
 
    
	// Fill interaction lists
    for (int i = 0; i < batches->num; i++) {

        pc_compute_interaction_list_remote2(0,
                tree_numpar, tree_radius,
                tree_x_mid, tree_y_mid, tree_z_mid,
                tree_num_children, tree_children,     
                batches_ind[i], batches_center[i], batches_radius[i],
                &(temp_tree_inter_list[i]), &(temp_direct_inter_list[i]),
                &(sizeof_tree_inter_list[i]), &(sizeof_direct_inter_list[i]),
                &(num_tree_inter[i]), &(num_direct_inter[i]));

        num_batch_approx[i] += num_tree_inter[i];
        num_batch_direct[i] += num_direct_inter[i];

    }


    for (int j = 0; j < batches->num; j++) {
        for (int i = 0; i < num_tree_inter[j]; i++) {

            int node_index = temp_tree_inter_list[j][i];
            approx_list_unpacked[node_index] = node_index;
        }

        for (int i = 0; i < num_direct_inter[j]; i++) {

            int node_index = temp_direct_inter_list[j][i];
            direct_list[node_index] = node_index;
        }
    }


    int approx_counter = 0;
    for (int i = 0; i < tree_numnodes; i++) {
        if (approx_list_unpacked[i] > -1) {
            approx_list_packed[approx_counter] = i;
            approx_counter++;
        }
    }


    free_matrix(temp_tree_inter_list);
    free_matrix(temp_direct_inter_list);

    free_vector(sizeof_tree_inter_list);
    free_vector(sizeof_direct_inter_list);
    
    free_vector(num_tree_inter);
    free_vector(num_direct_inter);

    return;

} /* END of function pc_treecode */


