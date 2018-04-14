#ifndef H_TREEFUNCTIONS_H
#define H_TREEFUNCTIONS_H

#include "tnode.h"


/* declaration of treecode support functions */

/* used by cluster-particle and particle-cluster */
void remove_node(struct tnode *p);

void cleanup(struct tnode *p);


/* used by cluster-particle and particle-cluster Coulomb */
void setup(double *x, double *y, double *z, int numpars,
           int order, double theta, double *xyzminmax);

void comp_tcoeff(double dx, double dy, double dz);


/* used by cluster-particle and particle-cluster Yukawa */
void setup_yuk(double *x, double *y, double *z, int numpars,
               int order, double theta, double *xyzminmax);

void comp_tcoeff_yuk(double dx, double dy, double dz, double kappa);


/* used by cluster-particle */
void cp_create_tree_n0(struct tnode **p, int ibeg, int iend,
                       double *x, double *y, double *z,
                       int maxparnode, double *xyzmm,
                       int level);

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
                 double *xS, double *yS, double *zS, double *qS,
                 double *xT, double *yT, double *zT, double *tpeng,
                 double *EnP, int numparsS, int numparsT,
                 double *timetree);

void compute_cp1(struct tnode *p, double *EnP,
                 double *x, double *y, double *z);

void cp_comp_direct(double *EnP, int ibeg, int iend,
                    double *x, double *y, double *z);


/* used by cluster-particle Yukawa */
void cp_treecode_yuk(struct tnode *p,
                     double *xS, double *yS, double *zS, double *qS,
                     double *xT, double *yT, double *zT, double *tpeng,
                     double *EnP, int numparsS, int numparsT,
                     double kappa, double *timetree);

void compute_cp1_yuk(struct tnode *p, double *EnP,
                     double *x, double *y, double *z,
                     double kappa);

void cp_comp_direct_yuk(double *EnP, int ibeg, int iend,
                        double *x, double *y, double *z,
                        double kappa);


/* used by particle-cluster */
void pc_create_tree_n0(struct tnode **p, int ibeg, int iend,
                       double *x, double *y, double *z, double *q,
                       int maxparnode, double *xyzmm,
                       int level);

void pc_partition_8(double *x, double *y, double *z, double *q,
                    double xyzmms[6][8], double xl, double yl, double zl,
                    double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid,
                    int ind[8][2]);

void pc_comp_ms(struct tnode *p, double *x, double *y, double *z, double *q);


/* used by particle-cluster Coulomb */
void pc_treecode(struct tnode *p, double *xS, double *yS, double *zS,
                 double *qS, double *xT, double *yT, double *zT,
                 double *tpeng, double *EnP, int numparsS, int numparsT,
                 int **batch_index, double **batch_center, double *batch_radius,
                 int batch_num, int *batch_reorder);

void compute_pc(struct tnode *p, double *EnP,
                double *x, double *y, double *z, double *q,
                double *xT, double *yT, double *zT,
                int *batch_ind, double *batch_mid, double batch_rad,
                int *batch_reorder);

void pc_comp_direct(double *EnP, int ibeg, int iend,
                    double *x, double *y, double *z, double *q,
                    int batch_ibeg, int batch_iend,
                    double *xT, double *yT, double *zT,
                    int *batch_reorder);


/* used by cluster-particle Yukawa */
void pc_treecode_yuk(struct tnode *p,
                     double *xS, double *yS, double *zS, double *qS,
                     double *xT, double *yT, double *zT, double *tpeng,
                     double *EnP, int numparsS, int numparsT,
                     double kappa);

void compute_pc_yuk(struct tnode *p, double *EnP,
                    double *x, double *y, double *z, double *q,
                    double kappa);

void pc_comp_direct_yuk(double *EnP, int ibeg, int iend,
                        double *x, double *y, double *z, double *q,
                        double kappa);


/* batch functions */
void cp_setup_batch(double *x, double *y, double *z,
           int numpars, int batch_size, double *xyzminmax, int **batch_reorder,
           int *batch_num, int ***batch_index, double ***batch_center,
           double **batch_radius);

void cp_create_batch(struct tnode **p, int ibeg, int iend,
                     double *x, double *y, double *z,
                     int maxparnode, double *xyzmm, int level,
                     int *batch_reorder, int *batch_num,
                     int **batch_index, double **batch_center,
                     double *batch_radius);

void cp_partition_batch(double *x, double *y, double *z, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *batch_reorder);
#endif /* H_TREEFUNCTIONS_H */
