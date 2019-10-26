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
					int *approx_list_unpacked,int *approx_list_packed, int *direct_list, int numnodes)
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

    for (i = 0; i < numnodes; i++) approx_list_unpacked[i] = -1;
    for (i = 0; i < numnodes; i++) approx_list_packed[i] = -1;
   	for (i = 0; i < numnodes; i++) direct_list[i] = -1;


    // Make interaction lists and set to -1
    int *temp_tree_inter_list, *temp_direct_inter_list;
    make_vector(temp_tree_inter_list, batches->num * numnodes);
	make_vector(temp_direct_inter_list, batches->num * numnodes);

    //size_t ii;
    //size_t loopsize = (size_t)batches->num * (size_t)numnodes;
    int loopsize = batches->num * numnodes;
	for (i = 0; i < loopsize; i++) temp_tree_inter_list[i] = -1;
	for (i = 0; i < loopsize; i++) temp_direct_inter_list[i] = -1;
 
    
	// Fill interaction lists
    for (i = 0; i < batches->num; i++) {
/*
        pc_compute_interaction_list_remote(tree_numnodes, tree_level, tree_numpar,
                tree_radius, tree_x_mid, tree_y_mid, tree_z_mid,
                batches_ind[i], batches_center[i], batches_radius[i],
                &(temp_tree_inter_list[i*numnodes]), &(temp_direct_inter_list[i*numnodes]));
*/

        int tree_index_counter = 0, direct_index_counter = 0;


        pc_compute_interaction_list_remote2(0,
                tree_numpar, tree_radius,
                tree_x_mid, tree_y_mid, tree_z_mid,
                tree_num_children, tree_children,     
                batches_ind[i], batches_center[i], batches_radius[i],
                &(temp_tree_inter_list[i*numnodes]), &(temp_direct_inter_list[i*numnodes]),
                &tree_index_counter, &direct_index_counter);

    }

//    printf("temp_tree_inter_list\n\n");
//        for (int i=0;i<numnodes;i++){
//        	printf("%i\n", temp_tree_inter_list[i]);
//	}
//    if (rank==0){
//	printf("temp_direct_inter_list\n\n");
//		for (int i=0;i<numnodes;i++){
//			printf("%i\n", temp_direct_inter_list[i]);
//		}
//    }


    // Update masks using interaction lists (overkill, but okay for now)
    int approx_counter=0, direct_counter=0;
    for (i=0; i<numnodes; i++){
    	for (j=0;j<batches->num;j++){
    		if (temp_tree_inter_list[j*numnodes+i]!=-1){ // then at least one target batch accepted the MAC for the ith node

    			approx_list_unpacked[i]=i;
    			approx_list_packed[approx_counter]=i;
				approx_counter+=1;
				break;
    		}
    	}
		for (j=0;j<batches->num;j++){
    		if (temp_direct_inter_list[j*numnodes+i]!=-1){ // then at least one target batch interacts directly with the ith node
//    			if (i==0) printf("Batch %i is putting the root in the direct list.\n", j);
    			direct_list[i]=i;
    			direct_counter+=1;
    			break;
			}
    	}
    }



/*
   	for (j=0;j<batches->num;j++){
        for (i=0; i<numnodes; i++){
    		if (temp_tree_inter_list[j*numnodes+i]=-1) { // then at least one target batch accepted the MAC for the ith node
                break;
            } else {
    			approx_list_unpacked[temp_tree_inter_list[j*numnodes+i]] = temp_tree_inter_list[j*numnodes+i];
    			//approx_list_packed[approx_counter]=i;
    		}
    	}

        for (i=0; i<numnodes; i++) {
    		if (temp_direct_inter_list[j*numnodes+i]=-1) { // then at least one target batch interacts directly with the ith node
                break;
            } else {
    			direct_list[temp_direct_inter_list[j*numnodes+i]] = temp_direct_inter_list[j*numnodes+i];
			}
    	}
    }

    int approx_counter=0, direct_counter=0;
    for (i=0; i<numnodes; i++) {
        if (approx_list_unpacked[i]!=-1) {
            approx_list_packed[approx_counter]=i;
            approx_counter++;
        }
    }
  */
  
//    exit(0);

/*
    if (rank==0){
    for (int i=0;i<numnodes*batches->num; i++){
               printf("batch %d, approx_list[%i] = %i\n", i/numnodes, i, temp_tree_inter_list[i]);
       }
    for (int i=0;i<numnodes*batches->num; i++){
       printf("batch %d, direct_list[%i] = %i\n", i/numnodes, i, temp_direct_inter_list[i]);
       }
    }
*/
    free_vector(temp_tree_inter_list);
    free_vector(temp_direct_inter_list);


    // Example:
    // At end of this function, approx_list and direct_list look like [c0, c3, c5, -1, -1, -1 ... ]
    // indicating the 0th, 3rd, and 5th clusters in remote tree array are needed.

//    printf("Exiting remote_interaction_lists\n");
    return;

} /* END of function pc_treecode */


