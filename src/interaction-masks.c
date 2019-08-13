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


void interaction_masks(const struct tnode_array *tree_array, struct batch *batches,
                              int *approx_mask, int *direct_mask)
{
    /* local variables */
    int i, j;

    int **batches_ind;
    double **batches_center;
    double *batches_radius;

    int tree_numnodes;
    const int *tree_numpar, *tree_level;
    const double *tree_x_mid, *tree_y_mid, *tree_z_mid, *tree_radius;

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

    for (i = 0; i < numnodes; i++){
    	approx_mask[i] = 0;
    	direct_mask[i] = 0;
    }


    // Make interaction lists and set to -1
    double *tree_inter_list, *direct_inter_list;
    make_vector(tree_inter_list, batches->num * numnodes);
	make_vector(direct_inter_list, batches->num * numleaves);

	for (i = 0; i < batches->num * numnodes; i++)
	        tree_inter_list[i] = -1;

	for (i = 0; i < batches->num * numleaves; i++)
		direct_inter_list[i] = -1;


	// Fill interaction lists
    for (i = 0; i < batches->num; i++){
    	// fill interaction lists
        pc_compute_interaction_list(tree_numnodes, tree_level, tree_numpar,
                tree_radius, tree_x_mid, tree_y_mid, tree_z_mid,
                batches_ind[i], batches_center[i], batches_radius[i],
                &(tree_inter_list[i*numnodes]), &(direct_inter_list[i*numleaves]));
    }

    // Update masks using interaction lists (overkill, but okay for now)
    for (i=0; i<numnodes; i++){
    	for (j=0;j<batches->num;j++){
    		if (tree_inter_list[j*numnodes]!=-1){ // then at least one target batch accepted the MAC for the ith node
    			approx_mask[i]=1
    		}
    		if (direct_inter_list[j*numnodes]!=-1){ // then at least one target batch interacts directly with the ith node
    			direct_mask[i]=1
			}
    	}

    }

    free_vector(tree_inter_list);
    free_vector(direct_inter_list);

    return;

} /* END of function pc_treecode */


void allocate_sources(struct particles *sources, newlength){

	sources->num = newlength;
	realloc_vector(sources->x, newlength);
	realloc_vector(sources->y, newlength);
	realloc_vector(sources->z, newlength);
	realloc_vector(sources->q, newlength);
	realloc_vector(sources->w, newlength);

    return;
} /* END of function allocate_sources */

void reallocate_sources(struct particles *sources, length){

	sources->num = length;
	make_vector(sources->x, length);
	make_vector(sources->y, length);
	make_vector(sources->z, length);
	make_vector(sources->q, length);
	make_vector(sources->w, length);

    return;
} /* END of function reallocate_sources */


void allocate_cluster(struct particles *clusters, length){

    make_vector(clusters->x, length);
    make_vector(clusters->y, length);
    make_vector(clusters->z, length);
    make_vector(clusters->q, length);
    make_vector(clusters->w, length);  // will be used in singularity subtraction
    clusters->num=length;

    return;
} /* END of function allocate_cluster */

void allocate_tree_array(tnode_array *let_tree_array, length){

	let_tree_array->numnodes = length;
	make_vector(let_tree_array->ibeg, length);
	make_vector(let_tree_array->iend, length);
	make_vector(let_tree_array->numpar, length);
	make_vector(let_tree_array->x_mid, length);
	make_vector(let_tree_array->y_mid, length);
	make_vector(let_tree_array->z_mid, length);
	make_vector(let_tree_array->x_min, length);
	make_vector(let_tree_array->y_min, length);
	make_vector(let_tree_array->z_min, length);
	make_vector(let_tree_array->x_max, length);
	make_vector(let_tree_array->y_max, length);
	make_vector(let_tree_array->z_max, length);
	make_vector(let_tree_array->level, length);
	make_vector(let_tree_array->cluster_ind, length);
	make_vector(let_tree_array->radius, length);

    return;
} /* END of function allocate_tree_array */


void reallocate_tree_array(tnode_array *let_tree_array, newlength){

	let_tree_array->numnodes = newlength;
	realloc_vector(let_tree_array->ibeg, newlength);
	realloc_vector(let_tree_array->iend, newlength);
	realloc_vector(let_tree_array->numpar, newlength);
	realloc_vector(let_tree_array->x_mid, newlength);
	realloc_vector(let_tree_array->y_mid, newlength);
	realloc_vector(let_tree_array->z_mid, newlength);
	realloc_vector(let_tree_array->x_min, newlength);
	realloc_vector(let_tree_array->y_min, newlength);
	realloc_vector(let_tree_array->z_min, newlength);
	realloc_vector(let_tree_array->x_max, newlength);
	realloc_vector(let_tree_array->y_max, newlength);
	realloc_vector(let_tree_array->z_max, newlength);
	realloc_vector(let_tree_array->level, newlength);
	realloc_vector(let_tree_array->cluster_ind, newlength);
	realloc_vector(let_tree_array->radius, newlength);

    return;
} /* END of function allocate_tree_array */


