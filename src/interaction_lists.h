#ifndef H_INTERACTIONLISTS_H
#define H_INTERACTIONLISTS_H

void Interaction_MakeList(const struct tnode_array *tree_array, struct tnode_array *batches,
                          int *tree_inter_list, int *direct_inter_list, int approx_offset,
                          int direct_offset,
                          int interpolationOrder, double sizeCheckFactor);

void Interaction_MakeListRemote(const struct tnode_array *tree_array, struct tnode_array *batches,
                                int *approx_mask_packed, int *approx_mask_unpacked, int *direct_mask,
                                int *num_batch_approx, int *num_batch_direct,
                                int interpolationOrder, double sizeCheckFactor);
                                
void Interaction_CP_MakeListRemote(const struct tnode_array *tree_array, struct tnode_array *batches,
                                   int *direct_list, int interpolationOrder, double sizeCheckFactor)

#endif /* H_INTERACTIONLISTS_H */
