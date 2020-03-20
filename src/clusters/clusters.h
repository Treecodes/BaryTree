#ifndef H_CLUSTERFUNCTIONS_H
#define H_CLUSTERFUNCTIONS_H

#include "../utilities/enums.h"

#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"

#include "struct_clusters.h"

void Clusters_PC_Setup(struct Clusters **clusters, struct Particles *sources, int order,
                       struct Tree *tree_array,
                       APPROXIMATION approxName, SINGULARITY singularity);

void Clusters_CP_Setup(struct Clusters **clusters, int order, struct Tree *tree_array,
                       APPROXIMATION approxName, SINGULARITY singularity);

void Clusters_Alloc(struct Clusters **clusters_addr, int length,
                    APPROXIMATION approxName, SINGULARITY singularity);

void Clusters_Free(struct Clusters *clusters);

void Clusters_Free_Win(struct Clusters *clusters);

#endif /* H_CLUSTERFUNCTIONS_H */
