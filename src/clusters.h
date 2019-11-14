#ifndef H_CLUSTERS_H
#define H_CLUSTERS_H

#include "nodes_struct.h"
#include "particles_struct.h"

void Clusters_PC_Setup(struct particles *clusters, struct particles *sources, int order,
                       struct tnode_array *tree_array, char *approxName, char *singularityHandling);

#endif /* H_CLUSTERS_H */
