#ifndef H_TREE_LINKED_LIST_FUNCTIONS_H
#define H_TREE_LINKED_LIST_FUNCTIONS_H

#include "../particles/struct_particles.h"
#include "struct_tree_linked_list_node.h"


void TreeLinkedList_Targets_Construct(struct TreeLinkedListNode **p, struct TreeLinkedListNode *parent,
                struct Particles *targets,
                int ibeg, int iend, int maxparnode, double *xyzmm, int *numnodes, int *numleaves,
                int *min_leaf_size, int *max_leaf_size, int *max_depth, int current_level);

void TreeLinkedList_Sources_Construct(struct TreeLinkedListNode **p, struct TreeLinkedListNode *parent,
                struct Particles *sources,
                int ibeg, int iend, int maxparnode, double *xyzmm,
                int *numnodes, int *numleaves, int *min_leaf_size, int *max_leaf_size,
                int *max_depth, int current_level);

int TreeLinkedList_SetIndex(struct TreeLinkedListNode *p, int index);

void TreeLinkedList_Free(struct TreeLinkedListNode **p_addr);


#endif /* H_TREE_LINKED_LIST_FUNCTIONS */
