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

void comp_tcoeff(double dx, double dy, double dz);


/* used by cluster-particle and particle-cluster Yukawa */
void setup_yuk(struct particles *particles, int order, double theta,
               double *xyzminmax);

void comp_tcoeff_yuk(double dx, double dy, double dz, double kappa);


/* used by cluster-particle */
void cp_create_tree_n0(struct tnode **p, struct particles *targets,
                       int ibeg, int iend, int maxparnode,
                       double *xyzmm, int level);

void cp_partition_8(double *x, double *y, double *z, double xyzmms[6][8],
                    double xl, double yl, double zl,
                    double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid,
                    int ind[8][2]);

void cp_comp_ms(struct tnode *p);

void compute_cp2(struct tnode *ap, double *x, double *y, double *z,
                 double *EnP);


/* used by cluster-particle Coulomb */
void cp_treecode(struct tnode *p,
                 struct particles *sources, struct particles *targets,
                 double *tpeng, double *EnP, double *timetree);

void compute_cp1(struct tnode *p, double *EnP,
                 double *x, double *y, double *z);

void cp_comp_direct(double *EnP, int ibeg, int iend,
                    double *x, double *y, double *z);


/* used by cluster-particle Yukawa */
void cp_treecode_yuk(struct tnode *p,
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

void pc_partition_8(double *x, double *y, double *z, double *q,
                    double xyzmms[6][8], double xl, double yl, double zl,
                    double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid,
                    int ind[8][2]);

void pc_comp_ms(struct tnode *p, double *x, double *y, double *z, double *q);


/* used by particle-cluster Coulomb */
void pc_treecode(struct tnode *p, struct batch *batches,
                 struct particles *sources, struct particles *targets,
                 double *tpeng, double *EnP);

void compute_pc(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS,
                double *xT, double *yT, double *zT, double *EnP);

void pc_comp_direct(int ibeg, int iend, int batch_ibeg, int batch_iend,
                    double *xS, double *yS, double *zS, double *qS,
                    double *xT, double *yT, double *zT, double *EnP);


/* used by cluster-particle Yukawa */
void pc_treecode_yuk(struct tnode *p, struct batch *batches,
                     struct particles *sources, struct particles *targets,
                     double kappa, double *tpeng, double *EnP);

void compute_pc_yuk(struct tnode *p,
                    int *batch_ind, double *batch_mid, double batch_rad,
                    double *xS, double *yS, double *zS, double *qS,
                    double *xT, double *yT, double *zT,
                    double kappa, double *EnP);

void pc_comp_direct_yuk(int ibeg, int iend, int batch_ibeg, int batch_iend,
                        double *xS, double *yS, double *zS, double *qS,
                        double *xT, double *yT, double *zT,
                        double kappa, double *EnP);


/* batch functions */
void setup_batch(struct batch **batches, double *batch_lim,
                 struct particles *particles, int batch_size);

void cp_create_batch(struct batch *batches, struct particles *particles,
                     int ibeg, int iend, int maxparnode, double *xyzmm);

void cp_partition_batch(double *x, double *y, double *z, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *batch_reorder);

void reorder_energies(int *batch_reorder, int numpars, double *tEn);
#endif /* H_TREEFUNCTIONS_H */
