#ifndef H_REMOTEINTERACTIONLISTS_H
#define H_REMOTEINTERACTIONLISTS_H

void remote_interaction_lists(const struct tnode_array *tree_array, struct batch *batches,
                              int *approx_mask_packed, int *approx_mask_unpacked, int *direct_mask,
                              int *num_batch_approx, int *num_batch_direct);

#endif /* H_REMOTEINTERACTIONLISTS_H */
