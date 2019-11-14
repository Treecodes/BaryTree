#ifndef H_TREEFUNCTIONS_H
#define H_TREEFUNCTIONS_H

#include "nodes_struct.h"
#include "particles_struct.h"


/* declaration of treecode support functions */


/* used by cluster-particle and particle-cluster */
void remove_node(struct tnode *p);

void cleanup(struct tnode *p);


/* used by cluster-particle and particle-cluster Coulomb */
void setup(struct particles *particles1, struct particles *particles2, int order, double theta,
           double *xyzminmax);




/* used by cluster-particle */
void cp_create_tree_n0(struct tnode **p, struct particles *targets,
                       int ibeg, int iend, int maxparnode, double *xyzmm,
                       int level, int *numnodes, int * numleaves);

void compute_cp2(struct tnode *ap, double *x, double *y, double *z,
                 double *EnP);



/* used by particle-cluster */
void pc_create_tree_n0(struct tnode **p, struct particles *sources,
                       int ibeg, int iend, int maxparnode, double *xyzmm,
                       int level, int *numnodes, int *numleaves);

int pc_set_tree_index(struct tnode *p, int index);

void pc_create_tree_array(struct tnode *p, struct tnode_array *tree_array);




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

#endif /* H_TREEFUNCTIONS_H */
