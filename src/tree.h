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
                       int shrink, int maxparnode, double *xyzmm,
                       int level);

void cp_create_tree_lv(struct tnode **p, int ibeg, int iend,
                       double *x, double *y, double *z,
                       int shrink, int treelevel, double *xyzmm,
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
                       int shrink, int maxparnode, double *xyzmm,
                       int level);

void pc_create_tree_lv(struct tnode **p, int ibeg, int iend,
                       double *x, double *y, double *z, double *q,
                       int shrink, int treelevel, double *xyzmm,
                       int level);

void pc_partition_8(double *x, double *y, double *z, double *q,
                    double xyzmms[6][8], double xl, double yl, double zl,
                    double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid,
                    int ind[8][2]);

void pc_comp_ms(struct tnode *p, double *x, double *y, double *z, double *q);


/* used by particle-cluster Coulomb */
void pc_treecode(struct tnode *p,
                 double *xS, double *yS, double *zS, double *qS,
                 double *xT, double *yT, double *zT, double *tpeng,
                 double *EnP, int numparsS, int numparsT);

void compute_pc(struct tnode *p, double *EnP,
                double *x, double *y, double *z, double *q);

void pc_comp_direct(double *EnP, int ibeg, int iend,
                    double *x, double *y, double *z, double *q);


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


/* use for grid treecode */
void setup_grid(double *xyzminmax, int *xyzdim, int *xyzind, 
                int order, double theta);

void cp_create_tree_n0_grid(struct tnode **p, int maxparnode, 
                            double *xyzmm, int *xyzdim, int *xyzind, int level);

void cp_partition_8_grid(double xyzmms[6][8], int xyzdims[3][8], int xyzinds[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid);

void cp_treecode_grid(struct tnode *p, double *xS, double *yS, double *zS, double *qS, 
                      double *tpeng, double *EnP, int numparsS, int numparsT,
                      double *timetree);

void compute_cp1_grid(struct tnode *p, double *EnP);

void cp_comp_direct_grid(struct tnode *p, double *EnP);

void compute_cp2_grid(struct tnode *ap, double *EnP);

#endif /* H_TREEFUNCTIONS_H */
