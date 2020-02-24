#ifndef H_INTERACTIONLISTS_H
#define H_INTERACTIONLISTS_H

void InteractionList_Make(const struct tnode_array *tree_array, struct tnode_array *batches,
                          int *tree_inter_list, int *direct_inter_list, int approx_offset,
                          int direct_offset,
                          int interpolationOrder, double sizeCheckFactor);

void InteractionList_PC_MakeRemote(const struct tnode_array *tree_array, struct tnode_array *batches,
                                int *approx_mask_packed, int *approx_mask_unpacked, int *direct_mask,
                                int *num_batch_approx, int *num_batch_direct,
                                int interpolationOrder, double sizeCheckFactor);
                                
void InteractionList_CP_MakeRemote(const struct tnode_array *tree_array, struct tnode_array *batches,
                                   int *direct_list, int interpolationOrder, double sizeCheckFactor);

void InteractionList_CC_Make(const struct tnode_array *source_tree_array,
                             struct tnode_array *target_tree_array,
                             int ***approx_inter_list_addr, int ***direct_inter_list_addr,
                             int interpolationOrder, double sizeCheckFactor);
                             
void InteractionList_CC_MakeRemote(const struct tnode_array *source_tree_array,
                                   struct tnode_array *target_tree_array,
                                   int *approx_list_unpacked,int *approx_list_packed, int *direct_list,
                                   int interpolationOrder, double sizeCheckFactor);
#endif /* H_INTERACTIONLISTS_H */
