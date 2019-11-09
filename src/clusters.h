#ifndef H_CLUSTERS_H
#define H_CLUSTERS_H

#include "tnode.h"
#include "particles.h"


void Clusters_PC_SetupLagrange(struct particles *clusters, struct particles *sources, int order, 
                               struct tnode_array *tree_array, char *singularityHandling);

void Clusters_PC_SetupHermite(struct particles *clusters, struct particles *sources, int order, 
                              struct tnode_array *tree_array, char *singularityHandling);

#endif /* H_CLUSTERS_H */
