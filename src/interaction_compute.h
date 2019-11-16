#ifndef H_INTERACTIONCOMPUTE_H
#define H_INTERACTIONCOMPUTE_H

#include "struct_nodes.h"


void compute_cp2(struct tnode *ap, double *x, double *y, double *z,
                 double *EnP);


void Interaction_PC_Compute(struct tnode_array *tree_array, struct tnode_array *batches,
                            int *tree_inter_list, int *direct_inter_list,
                            double *xS, double *yS, double *zS, double *qS, double *wS,
                            double *xT, double *yT, double *zT, double *qT,
                            double *xC, double *yC, double *zC, double *qC, double *wC,
                            double *pointwisePotential, int interpolationOrder,
                            int numSources, int numTargets, int numClusters,
                            int offset_approx, int offset_direct,
                            char *kernelName, double kernel_parameter, char *singularityHandling,
                            char *approximationName);


void Interaction_Direct_Compute(double *source_x, double *source_y, double *source_z,
                            double *source_q, double *source_w,
                            double *target_x, double *target_y, double *target_z, double *target_q,
                            double *totalPotential, int numSources, int numTargets,
                            char *kernelName, double kernel_parameter, char *singularityHandling,
                            char *approximationName);


void Interaction_SubtractionPotentialCorrection(double *pointwisePotential, double *target_q, int numTargets,
                            char *kernelName, double kernel_parameter, char *singularityHandling);

#endif /* H_INTERACTIONCOMPUTE_H */
