#ifndef H_INTERACTION_COMPUTE_H
#define H_INTERACTION_COMPUTE_H

#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"
#include "../clusters/struct_clusters.h"
#include "../run_params/struct_run_params.h"
#include "../interaction_lists/struct_interaction_lists.h"


void InteractionCompute_PC(double *potential, struct Tree *tree, struct Tree *batches,
                           struct InteractionLists *interaction_list,
                           struct Particles *sources, struct Particles *targets,
                           struct Clusters *clusters, struct RunParams *run_params);


void InteractionCompute_CP(double *potential, struct Tree *tree, struct Tree *batches,
                           struct InteractionLists *interaction_list,
                           struct Particles *sources, struct Particles *targets,
                           struct Clusters *clusters, struct RunParams *run_params);


void InteractionCompute_CC(double *potential, struct Tree *source_tree, struct Tree *target_tree,
                           struct InteractionLists *interaction_list,
                           struct Particles *sources, struct Particles *targets,
                           struct Clusters *source_clusters, struct Clusters *target_clusters,
                           struct RunParams *run_params);


void InteractionCompute_Downpass(double *potential, struct Tree *tree,
                           struct Particles *targets, struct Clusters *clusters,
                           struct RunParams *run_params);


void InteractionCompute_Direct(double *potential,
                           struct Particles *sources, struct Particles *targets,
                           struct RunParams *run_params);


void InteractionCompute_SubtractionPotentialCorrection(double *potential, 
                           struct Particles *targets, struct RunParams *run_params);


#endif /* H_INTERACTION_COMPUTE_H */
