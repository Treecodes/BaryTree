#ifndef H_INTERACTIONCOMPUTE_H
#define H_INTERACTIONCOMPUTE_H

#include "struct_nodes.h"
#include "struct_kernel.h"


void InteractionCompute_PC(struct tnode_array *tree_array, struct tnode_array *batches,
                            int *tree_inter_list, int *direct_inter_list,
                            double *xS, double *yS, double *zS, double *qS, double *wS,
                            double *xT, double *yT, double *zT, double *qT,
                            double *xC, double *yC, double *zC, double *qC, double *wC,
                            double *pointwisePotential, int interpolationOrder,
                            int numSources, int numTargets, int numClusters,
                            int offset_approx, int offset_direct,
                            struct kernel *kernel, char *singularityHandling,
                            char *approximationName);


void InteractionCompute_CP_1(struct tnode_array *tree_array, struct tnode_array *batches,
                            int *tree_inter_list, int *direct_inter_list,
                            double *source_x, double *source_y, double *source_z,
                            double *source_charge, double *source_weight,
                            double *target_x, double *target_y, double *target_z, double *target_charge,
                            double *cluster_x, double *cluster_y, double *cluster_z,
                            double *cluster_charge, double *cluster_weight,
                            double *pointwisePotential, int interpolationOrder,
                            int numSources, int numTargets, int totalNumberOfInterpolationPoints,
                            int batch_approx_offset, int batch_direct_offset,
                            struct kernel *kernel, char *singularityHandling,
                            char *approximationName);
                            
                            
void InteractionCompute_CP_2(struct tnode_array *tree_array,
                            double *target_x, double *target_y, double *target_z, double *target_charge,
                            double *cluster_x, double *cluster_y, double *cluster_z,
                            double *cluster_charge, double *cluster_weight,
                            double *pointwisePotential, int interpolationOrder,
                            int numTargets, int totalNumberInterpolationPoints,
                            int totalNumberInterpolationCharges, int totalNumberInterpolationWeights,
                            char *singularityHandling, char *approximationName);


void InteractionCompute_Direct(double *source_x, double *source_y, double *source_z,
                            double *source_q, double *source_w,
                            double *target_x, double *target_y, double *target_z, double *target_q,
                            double *totalPotential, int numSources, int numTargets,
                            struct kernel *kernel, char *singularityHandling,
                            char *approximationName);


void InteractionCompute_SubtractionPotentialCorrection(double *pointwisePotential, double *target_q, int numTargets,
                            struct kernel *kernel, char *singularityHandling);

#endif /* H_INTERACTIONCOMPUTE_H */
