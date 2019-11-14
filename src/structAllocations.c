#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "array.h"
#include "tools.h"
#include "struct_particles.h"
#include "struct_nodes.h"



void allocate_sources(struct particles *sources, int length)  {

	sources->num = length;
	make_vector(sources->x, length);
	make_vector(sources->y, length);
	make_vector(sources->z, length);
	make_vector(sources->q, length);
	make_vector(sources->w, length);

    return;
}   /* END of function allocate_sources */


void allocate_cluster(struct particles *clusters, int length)  {

    make_vector(clusters->x, length);
    make_vector(clusters->y, length);
    make_vector(clusters->z, length);
    make_vector(clusters->q, length);
    make_vector(clusters->w, length);  // will be used in singularity subtraction
    clusters->num = length;

    return;
}   /* END of function allocate_cluster */


void allocate_tree_array(struct tnode_array *tree_array, int length)  {

	tree_array->numnodes = length;
	make_vector(tree_array->ibeg, length);
	make_vector(tree_array->iend, length);
	make_vector(tree_array->numpar, length);
	make_vector(tree_array->x_mid, length);
	make_vector(tree_array->y_mid, length);
	make_vector(tree_array->z_mid, length);
	make_vector(tree_array->x_min, length);
	make_vector(tree_array->y_min, length);
	make_vector(tree_array->z_min, length);
	make_vector(tree_array->x_max, length);
	make_vector(tree_array->y_max, length);
	make_vector(tree_array->z_max, length);
	make_vector(tree_array->level, length);
	make_vector(tree_array->cluster_ind, length);
	make_vector(tree_array->radius, length);

	make_vector(tree_array->num_children, length);
	make_vector(tree_array->children, 8*length);

    return;
}   /* END of function allocate_tree_array */

void free_tree_array(struct tnode_array *tree_array)  {

	free_vector(tree_array->ibeg);
	free_vector(tree_array->iend);
	free_vector(tree_array->numpar);
	free_vector(tree_array->x_mid);
	free_vector(tree_array->y_mid);
	free_vector(tree_array->z_mid);
	free_vector(tree_array->x_min);
	free_vector(tree_array->y_min);
	free_vector(tree_array->z_min);
	free_vector(tree_array->x_max);
	free_vector(tree_array->y_max);
	free_vector(tree_array->z_max);
	free_vector(tree_array->level);
	free_vector(tree_array->cluster_ind);
	free_vector(tree_array->radius);

	free_vector(tree_array->num_children);
	free_vector(tree_array->children);
    
    free(tree_array);
	tree_array = NULL;

    return;
}   /* END of function allocate_tree_array */


void reallocate_tree_array(struct tnode_array *tree_array, int newlength)  {

	tree_array->numnodes = newlength;
	realloc_vector(tree_array->ibeg, newlength);
	realloc_vector(tree_array->iend, newlength);
	realloc_vector(tree_array->numpar, newlength);
	realloc_vector(tree_array->x_mid, newlength);
	realloc_vector(tree_array->y_mid, newlength);
	realloc_vector(tree_array->z_mid, newlength);
	realloc_vector(tree_array->x_min, newlength);
	realloc_vector(tree_array->y_min, newlength);
	realloc_vector(tree_array->z_min, newlength);
	realloc_vector(tree_array->x_max, newlength);
	realloc_vector(tree_array->y_max, newlength);
	realloc_vector(tree_array->z_max, newlength);
	realloc_vector(tree_array->level, newlength);
	realloc_vector(tree_array->cluster_ind, newlength);
	realloc_vector(tree_array->radius, newlength);

	realloc_vector(tree_array->num_children, newlength);
	realloc_vector(tree_array->children, 8*newlength);

	return;

}   /* END of function allocate_tree_array */

