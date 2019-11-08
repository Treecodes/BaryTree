#ifndef H_TREEFUNCTIONS_H
#define H_TREEFUNCTIONS_H

#include "tnode.h"
#include "batch.h"
#include "particles.h"


/* declaration of treecode support functions */


/* used by cluster-particle and particle-cluster */
void remove_node(struct tnode *p);

void cleanup(struct tnode *p);


/* used by cluster-particle and particle-cluster Coulomb */
void setup(struct particles *particles, int order, double theta,
           double *xyzminmax);

void addNodeToArray_hermite(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster);
void addNodeToArray_hermite_SS(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster);
//void addNodeToArray_SS(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster);



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




void pc_interaction_list_treecode(struct tnode_array *tree_array, struct batch *batches,
								  int *tree_inter_list, int *direct_inter_list,
								  double *xS, double *yS, double *zS, double *qS, double *wS,
								  double *xT, double *yT, double *zT, double *qT,
								  double *xC, double *yC, double *zC, double *qC, double *wC,
								  double *totalPotential, double *pointwisePotential, int interpolationOrder,
								  int numSources, int numTargets, int numClusters,
                                  int offset_approx, int offset_direct,
								  char *kernelName, double kernel_parameter, char *singularityHandling,
								  char *approximationName);


/* used by particle-cluster Coulomb */

void pc_treecode_hermite_coulomb_SS(struct tnode *p, struct batch *batches,
                 struct particles *sources, struct particles *targets, struct particles *clusters,
				 double kappa, double *tpeng, double *EnP);

void compute_pc_hermite_SS(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double kappaSq, double *EnP,
				double * clusterX, double * clusterY, double * clusterZ, double * clusterQ,
				double * clusterQx,double * clusterQy,double * clusterQz,double * clusterQxy,double * clusterQyz,double * clusterQxz,double * clusterQxyz,
				double * clusterW, double * clusterWx,double * clusterWy,double * clusterWz,double * clusterWxy,double * clusterWyz,double * clusterWxz,double * clusterWxyz);



/* batch functions */
void setup_batch(struct batch **batches, double *batch_lim,
                 struct particles *particles, int batch_size);

void create_target_batch(struct batch *batches, struct particles *particles,
                     int ibeg, int iend, int maxparnode, double *xyzmm);

void cp_partition_batch(double *x, double *y, double *z, double *q, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *batch_reorder);

void create_source_batch(struct batch *batches, struct particles *particles,
                     int ibeg, int iend, int maxparnode, double *xyzmm);

void pc_partition_batch(double *x, double *y, double *z, double *q, double *w, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *batch_reorder);

void reorder_energies(int *batch_reorder, int numpars, double *tEn);
#endif /* H_TREEFUNCTIONS_H */
