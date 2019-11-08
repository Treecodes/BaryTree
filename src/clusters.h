#ifndef H_CLUSTERS_H
#define H_CLUSTERS_H

#include "tnode.h"
#include "batch.h"
#include "particles.h"


void Clusters_PC_SetupLagrange(struct particles *clusters, struct particles *sources,
                               struct tnode *troot, int order, 
                               struct tnode_array *tree_array, char *singularityHandling);

void fill_in_cluster_data_hermite(struct particles *clusters, struct particles *sources, struct tnode *troot, int order);

#endif /* H_CLUSTERS_H */
