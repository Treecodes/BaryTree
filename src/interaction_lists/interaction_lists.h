#ifndef H_INTERACTION_LISTS_H
#define H_INTERACTION_LISTS_H

#include "../run_params/struct_run_params.h"
#include "struct_interaction_lists.h"


void InteractionList_Make(struct InteractionList **interaction_list_addr,
                          const struct tnode_array *tree_array, struct tnode_array *batches,
                          struct RunParams *run_params);
                          
void InteractionList_Free(struct InteractionList *interaction_list_addr);

void InteractionList_PC_MakeRemote(const struct tnode_array *tree_array, struct tnode_array *batches,
                          int *approx_list_packed, int *approx_list_unpacked, int *direct_list,
                          struct RunParams *run_params);
                                
void InteractionList_CP_MakeRemote(const struct tnode_array *tree_array, struct tnode_array *batches,
                          int *direct_list, struct RunParams *run_params);


void InteractionList_CC_Make(const struct tnode_array *source_tree_array,
                          struct tnode_array *target_tree_array,
                          int ***approx_inter_list_addr, int ***direct_inter_list_addr,
                          struct RunParams *run_params);
                             
void InteractionList_CC_MakeRemote(const struct tnode_array *source_tree_array,
                          struct tnode_array *target_tree_array,
                          int *approx_list_unpacked,int *approx_list_packed, int *direct_list,
                          struct RunParams *run_params);


#endif /* H_INTERACTION_LISTS_H */
