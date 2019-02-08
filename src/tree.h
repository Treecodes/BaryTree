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

void fill_in_cluster_data(struct particles *clusters, struct particles *sources, struct tnode *troot, int order);
void fill_in_cluster_data_SS(struct particles *clusters, struct particles *sources, struct tnode *troot, int order);

void addNodeToArray(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster);
void addNodeToArray_SS(struct tnode *p, struct particles *sources, struct particles *clusters, int order, int numInterpPoints, int pointsPerCluster);
void comp_tcoeff(double dx, double dy, double dz);


/* used by cluster-particle and particle-cluster Yukawa */
void setup_yuk(struct particles *particles, int order, double theta,
               double *xyzminmax);

void comp_tcoeff_yuk(double dx, double dy, double dz, double kappa);


/* used by cluster-particle */
void cp_create_tree_n0(struct tnode **p, struct particles *targets,
                       int ibeg, int iend, int maxparnode, double *xyzmm,
                       int level);

void cp_partition_8(double *x, double *y, double *z, double *q, double xyzmms[6][8],
                    double xl, double yl, double zl,
                    double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid,
                    int ind[8][2]);

void cp_comp_ms(struct tnode *p);

void compute_cp2(struct tnode *ap, double *x, double *y, double *z,
                 double *EnP);


/* used by cluster-particle Coulomb */
void cp_treecode(struct tnode *p, struct batch *batches,
                 struct particles *sources, struct particles *targets,
                 double *tpeng, double *EnP, double *timetree);

void compute_cp1(struct tnode *p, double *EnP,
                 double *x, double *y, double *z);

void cp_comp_direct(double *EnP, int ibeg, int iend,
                    double *x, double *y, double *z);


/* used by cluster-particle Yukawa */
void cp_treecode_yuk(struct tnode *p, struct batch *batches,
                     struct particles *sources, struct particles *targets,
                     double kappa, double *tpeng, double *EnP,
                     double *timetree);


void compute_cp1_yuk(struct tnode *p, double *EnP,
                     double *x, double *y, double *z,
                     double kappa);

void cp_comp_direct_yuk(double *EnP, int ibeg, int iend,
                        double *x, double *y, double *z,
                        double kappa);


/* used by particle-cluster */
void pc_create_tree_n0(struct tnode **p, struct particles *sources,
                       int ibeg, int iend, int maxparnode, double *xyzmm,
                       int level);

void pc_create_tree_array(struct tnode *p, struct tnode_array *tree_array);

void pc_partition_8(double *x, double *y, double *z, double *q, double *w,
                    double xyzmms[6][8], double xl, double yl, double zl,
                    double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid,
                    int ind[8][2]);

void pc_comp_ms(struct tnode *p, double *x, double *y, double *z, double *q, double *w, double *clusterQ);
void pc_comp_ms_SS(struct tnode *p, double *x, double *y, double *z, double *q, double *w, double *clusterQ, double *clusterQ2);

void pc_comp_ms_gpu(struct tnode *p, double __restrict__ *xS, double __restrict__ *yS, double __restrict__ *zS, double __restrict__ *qS, double __restrict__ *wS,
		double __restrict__ *clusterX, double __restrict__ *clusterY, double __restrict__ *clusterZ, double __restrict__ *clusterQ);

void pc_comp_weights(struct tnode *p);



void pc_make_interaction_list(struct tnode *p, struct batch *batches,
                              int **tree_inter_list, int **direct_inter_list);

void pc_compute_interaction_list(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                int *batch_tree_list, int *batch_direct_list,
                int *tree_index_counter, int *direct_index_counter);


/* used by particle-cluster Coulomb */
void pc_treecode(struct tnode *p, struct batch *batches,
                 struct particles *sources, struct particles *targets, struct particles *clusters,
                 double *tpeng, double *EnP);

void compute_pc(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double *EnP,
				double *clusterX, double *clusterY, double *clusterZ, double *clusterM);

void pc_comp_direct(int ibeg, int iend, int batch_ibeg, int batch_iend,
                    double *xS, double *yS, double *zS, double *qS, double *wS,
                    double *xT, double *yT, double *zT, double *qT, double *EnP);


/* used by particle-cluster Yukawa */
void pc_treecode_yuk(struct tnode *p, struct batch *batches,
                     struct particles *sources, struct particles *targets, struct particles *clusters,
                     double kappa, double *tpeng, double *EnP);

void compute_pc_yuk(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double kappa, double *EnP,
				double * clusterX, double * clusterY, double * clusterZ, double * clusterM );

void pc_comp_direct_yuk(int ibeg, int iend, int batch_ibeg, int batch_iend,
                        double *xS, double *yS, double *zS, double *qS, double *wS,
                        double *xT, double *yT, double *zT, double *qT,
                        double kappa, double *EnP);

/* used by particle-cluster Yukawa w/ singularity subtraction */

void pc_treecode_yuk_SS(struct tnode *p, struct batch *batches,
                     struct particles *sources, struct particles *targets, struct particles *clusters,
                     double kappa, double *tpeng, double *EnP);

void compute_pc_yuk_SS(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double kappa, double *EnP,
				double * clusterX, double * clusterY, double * clusterZ, double * clusterM , double * clusterM2);

void pc_comp_direct_yuk_SS(int ibeg, int iend, int batch_ibeg, int batch_iend,
                    double *xS, double *yS, double *zS, double *qS, double *wS,
                    double *xT, double *yT, double *zT, double *qT, double kappa, double *EnP);


/* used by particle-cluster Coulomb kernel w/ singularity subtraction */

void pc_treecode_coulomb_SS(struct tnode *p, struct batch *batches,
                     struct particles *sources, struct particles *targets, struct particles *clusters,
                     double kappaSq, double *tpeng, double *EnP);

void compute_pc_coulomb_SS(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double kappaSq, double *EnP,
				double * clusterX, double * clusterY, double * clusterZ, double * clusterM, double * clusterM2 );

void pc_comp_direct_coulomb_SS(int ibeg, int iend, int batch_ibeg, int batch_iend,
                    double *xS, double *yS, double *zS, double *qS, double *wS,
                    double *xT, double *yT, double *zT, double *qT, double kappaSq, double *EnP);


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
