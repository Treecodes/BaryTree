#ifndef H_TREEFUNCTIONS_H
#define H_TREEFUNCTIONS_H

#include "tnode.h"

/* declaration of treecode support functions*/

void remove_node(struct tnode *p);

void cleanup(struct tnode *p);

void comp_direct(double *EnP, int ibeg, int iend,
                 double *x, double *y, double *z,
                 int arrdim);

void comp_direct_yuk(double *EnP, int ibeg, int iend,
                     double *x, double *y, double *z,
                     int arrdim, double kappa);

void comp_cms(struct tnode *p);

void comp_tcoeff(double dx, double dy, double dz);

void comp_tcoeff_yuk(double dx, double dy, double dz,
                     double kappa);

void compute_cp2(struct tnode *ap, 
                 double *x, double *y, double *z,
                 double *EnP, int arrdim);

void compute_cp1(struct tnode *p, double *EnP,
                 double *x, double *y, double *z,
                 int arrdim);

void compute_cp1_yuk(struct tnode *p, double *EnP,
                     double *x, double *y, double *z,
                     int arrdim, double kappa);

void cp_treecode(struct tnode *p,
                 double *xS, double *yS, double *zS, double *qS,
                 double *xT, double *yT, double *zT, double *tpeng,
                 double *EnP, int numparsS, int numparsT);

void cp_treecode_yuk(struct tnode *p,
                     double *xS, double *yS, double *zS, double *qS,
                     double *xT, double *yT, double *zT, double *tpeng,
                     double *EnP, int numparsS, int numparsT,
                     double kappa);

void partition_8(double *x, double *y, double *z, double **xyzmms,
                 double xl, double yl, double zl, 
                 double lmax, int *numposchild,
                 double x_mid, double y_mid, double z_mid,
                 int **ind, int arrdim);

void create_tree_lv(struct tnode **p, int ibeg, int iend,
                    double *x, double *y, double *z,
                    int shrink, int treelevel, double *xyzmm,
                    int level, int arrdim);

void create_tree_n0(struct tnode **p, int ibeg, int iend,
                    double *x, double *y, double *z,
                    int shrink, int maxparnode, double *xyzmm,
                    int level, int arrdim);

void setup(double *x, double *y, double *z, int numpars,
           int order, double theta, double *xyzminmax);

void setup_yuk(double *x, double *y, double *z, int numpars,
               int order, double theta, double *xyzminmax);

#endif /* H_TREEFUNCTIONS_H */
