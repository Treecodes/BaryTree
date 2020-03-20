#ifndef H_INTERACTIONCOMPUTE_H
#define H_INTERACTIONCOMPUTE_H

#include "../tree/struct_tree.h"
#include "../run_params/struct_run_params.h"
#include "../interaction_lists/struct_interaction_lists.h"


void InteractionCompute_PC(struct Tree *tree, struct Tree *batches,
                           struct InteractionLists *interaction_list,
                           double *xS, double *yS, double *zS, double *qS, double *wS,
                           double *xT, double *yT, double *zT, double *qT,
                           double *xC, double *yC, double *zC, double *qC, double *wC,
                           double *pointwisePotential,
                           int numSources, int numTargets, int totalNumberOfInterpolationPoints,
                           struct RunParams *run_params);


void InteractionCompute_CP(struct Tree *tree, struct Tree *batches,
                           struct InteractionLists *interaction_list,
                           double *source_x, double *source_y, double *source_z,
                           double *source_charge, double *source_weight,
                           double *target_x, double *target_y, double *target_z, double *target_charge,
                           double *cluster_x, double *cluster_y, double *cluster_z,
                           double *cluster_charge, double *cluster_weight,
                           double *pointwisePotential,
                           int numSources, int numTargets, int totalNumberOfInterpolationPoints,
                           struct RunParams *run_params);
                          
                          
void InteractionCompute_CC(struct Tree *source_tree, struct Tree *target_tree,
                           struct InteractionLists *interaction_list,
                           double *source_x, double *source_y, double *source_z,
                           double *source_q, double *source_w,
                           double *target_x, double *target_y, double *target_z, double *target_q,
                           double *source_cluster_x, double *source_cluster_y, double *source_cluster_z,
                           double *source_cluster_q, double *source_cluster_w,
                           double *target_cluster_x, double *target_cluster_y, double *target_cluster_z,
                           double *target_cluster_q, double *target_cluster_w,
                           double *pointwisePotential,
                           int numSources, int numTargets, int numSourceClusterPoints, int numTargetClusterPoints,
                           struct RunParams *run_params);


void InteractionCompute_Downpass(struct Tree *tree,
                           double *target_x, double *target_y, double *target_z, double *target_charge,
                           double *cluster_x, double *cluster_y, double *cluster_z,
                           double *cluster_charge, double *cluster_weight,
                           double *pointwisePotential, int numTargets,
                           int totalNumberInterpolationCharges, int totalNumberInterpolationWeights,
                           struct RunParams *run_params);


void InteractionCompute_Direct(double *source_x, double *source_y, double *source_z,
                           double *source_q, double *source_w,
                           double *target_x, double *target_y, double *target_z, double *target_q,
                           double *totalPotential, int numSources, int numTargets,
                           struct RunParams *run_params);


void InteractionCompute_SubtractionPotentialCorrection(double *pointwisePotential, double *target_q, int numTargets,
                           struct RunParams *run_params);


#endif /* H_INTERACTIONCOMPUTE_H */
