#ifndef H_INTERACTION_LISTS_H
#define H_INTERACTION_LISTS_H

#include "../tree/struct_tree.h"
#include "../run_params/struct_run_params.h"
#include "struct_interaction_lists.h"


void InteractionLists_Make(struct InteractionLists **interaction_list_addr,
                          const struct Tree *source_tree, const struct Tree *target_tree,
                          const struct RunParams *run_params);
                          
void InteractionLists_Free(struct InteractionLists **interaction_list_addr);

void InteractionLists_MakeRemote(const struct Tree *source_tree, const struct Tree *target_tree,
                          int *approx_list_packed, int *approx_list_unpacked, int *direct_list,
                          const struct RunParams *run_params);


#endif /* H_INTERACTION_LISTS_H */
