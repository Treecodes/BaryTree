#ifndef H_INTERACTION_LISTS_H
#define H_INTERACTION_LISTS_H

#include "../run_params/struct_run_params.h"
#include "struct_interaction_lists.h"


void InteractionLists_Make(struct InteractionLists **interaction_list_addr,
                          const struct tnode_array *source_tree_array,
                          const struct tnode_array *target_tree_array,
                          struct RunParams *run_params);
                          
void InteractionLists_Free(struct InteractionLists *interaction_list_addr);


void InteractionLists_MakeRemote(const struct tnode_array *source_tree_array,
                          const struct tnode_array *target_tree_array,
                          int *approx_list_packed, int *approx_list_unpacked, int *direct_list,
                          struct RunParams *run_params);


#endif /* H_INTERACTION_LISTS_H */
