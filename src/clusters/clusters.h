#ifndef H_CLUSTER_FUNCTIONS_H
#define H_CLUSTER_FUNCTIONS_H

#include "../utilities/enums.h"

#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"

#include "struct_clusters.h"


void Clusters_Sources_Construct(struct Clusters **clusters, const struct Particles *sources,
                const struct Tree *tree, const struct RunParams *run_params);

void Clusters_Targets_Construct(struct Clusters **clusters, const struct Particles *targets,
                const struct Tree *tree, const struct RunParams *run_params);

void Clusters_Alloc(struct Clusters **clusters_addr, int length,
                const struct RunParams *run_params);

void Clusters_Free(struct Clusters **clusters_addr);

void Clusters_Free_Win(struct Clusters **clusters_addr);


#endif /* H_CLUSTER_FUNCTIONS_H */
