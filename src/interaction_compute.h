#ifndef H_INTERACTIONCOMPUTE_H
#define H_INTERACTIONCOMPUTE_H

#include "struct_nodes.h"


void compute_cp2(struct tnode *ap, double *x, double *y, double *z,
                 double *EnP);


void pc_interaction_list_treecode(struct tnode_array *tree_array, struct tnode_array *batches,
								  int *tree_inter_list, int *direct_inter_list,
								  double *xS, double *yS, double *zS, double *qS, double *wS,
								  double *xT, double *yT, double *zT, double *qT,
								  double *xC, double *yC, double *zC, double *qC, double *wC,
								  double *totalPotential, double *pointwisePotential, int interpolationOrder,
								  int numSources, int numTargets, int numClusters,
                                  int offset_approx, int offset_direct,
								  char *kernelName, double kernel_parameter, char *singularityHandling,
								  char *approximationName);

#endif /* H_INTERACTIONCOMPUTE_H */
