#ifndef H_TREEDRIVER_H
#define H_TREEDRIVER_H

void interaction_masks(const struct tnode_array *tree_array, struct batch *batches,
                              int *approx_mask, int *direct_mask);
