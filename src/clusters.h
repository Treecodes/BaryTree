#ifndef H_CLUSTERFUNCTIONS_H
#define H_CLUSTERFUNCTIONS_H

#include "struct_nodes.h"
#include "struct_particles.h"
#include "struct_clusters.h"

void Clusters_PC_Setup(struct clusters **clusters, struct particles *sources, int order,
                       struct tnode_array *tree_array, char *approxName, char *singularityHandling);

void Clusters_Alloc(struct clusters *clusters, int length, char *approxName, char *singularityHandling);

void Clusters_Free(struct clusters *clusters);

void Clusters_Free_Win(struct clusters *clusters);

#endif /* H_CLUSTERFUNCTIONS_H */
