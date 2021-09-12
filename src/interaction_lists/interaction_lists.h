#ifndef H_INTERACTION_LISTS_H
#define H_INTERACTION_LISTS_H

#include "../tree/struct_tree.h"
#include "../run_params/struct_run_params.h"
#include "struct_interaction_lists.h"


void InteractionLists_Make(struct InteractionLists **interaction_list_addr,
                          const struct Tree *source_tree, const struct Tree *target_tree,
                          const struct RunParams *run_params);
                          
void InteractionLists_Free(struct InteractionLists **interaction_list_addr);


#endif /* H_INTERACTION_LISTS_H */
