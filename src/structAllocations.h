#ifndef H_STRUCTALLOCATIONS_H
#define H_STRUCTALLOCATIONS_H

#include "array.h"
#include "treedriver.h"
#include "tools.h"
#include "particles.h"
#include "sort.h"
#include "tnode.h"
#include "batch.h"

/*
 * declaration of struct allocation functions
 *
 * partition determines the index MIDIND, after partitioning in place the arrays a, b, c,
 * and q, such that a(ibeg:midind) <= val and a(midind+1:iend) > val. If on entry, ibeg >
 * iend, or a(ibeg:iend) > val then midind is returned as ibeg-1.
 */

void reallocate_sources(struct particles *sources, int newlength);

void allocate_sources(struct particles *sources, int length);


void allocate_cluster(struct particles *clusters, int length);


void reallocate_cluster(struct particles *clusters, int newlength);

void allocate_tree_array(struct tnode_array *let_tree_array, int length);

void free_tree_array(struct tnode_array *let_tree_array);


void reallocate_tree_array(struct tnode_array *let_tree_array, int newlength);   /* END of function allocate_tree_array */


#endif /* H_STRUCTALLOCATIONS_H */
