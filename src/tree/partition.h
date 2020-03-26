#ifndef H_PARTITION_H
#define H_PARTITION_H

/* 
 * declaration of partition functions
 *
 * partition determines the index MIDIND, after partitioning in place the arrays a, b, c,
 * and q, such that a(ibeg:midind) <= val and a(midind+1:iend) > val. If on entry, ibeg >
 * iend, or a(ibeg:iend) > val then midind is returned as ibeg-1.
 */

void cp_partition(double *a, double *b, double *c, double *d, int *indarr,
                  int ibeg, int iend, double val, int *midind);

void pc_partition(double *a, double *b, double *c, double *d, double *w, int *indarr,
                  int ibeg, int iend, double val, int *midind);

void pc_partition_8(double *x, double *y, double *z, double *q, double *w, int *indarr,
                    double xyzmms[6][8], double xl, double yl, double zl,
                    int *numposchild, int max_num_children,
                    double x_mid, double y_mid, double z_mid,
                    int ind[8][2]);

void cp_partition_8(double *x, double *y, double *z, double *q, int *indarr,
                    double xyzmms[6][8], double xl, double yl, double zl,
                    int *numposchild, int max_num_children,
                    double x_mid, double y_mid, double z_mid,
                    int ind[8][2]);


#endif /* H_PARTITION_H */
