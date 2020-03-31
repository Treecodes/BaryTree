#include <stdio.h>
#include <math.h>

#include "../utilities/tools.h"

#include "partition.h"

/* 
 * partition determines the index MIDIND, after partitioning in place the arrays a, b, c,
 * and q, such that a(ibeg:midind) <= val and a(midind+1:iend) > val. If on entry, ibeg >
 * iend, or a(ibeg:iend) > val then midind is returned as ibeg-1.
 */

void cp_partition(double *a, double *b, double *c, double *d, int *indarr,
                  int ibeg, int iend, double val, int *midind)
{
        double ta, tb, tc, td;
        int lower, upper, tind;

        if (ibeg < iend) {
        /* 
         * temporarily stores ibeg entries and set a(ibeg) = val
         * for the partitioning algorithm.
         */
         
                ta = a[ibeg - 1];
                tb = b[ibeg - 1];
                tc = c[ibeg - 1];
                td = d[ibeg - 1];
                tind = indarr[ibeg - 1];
                a[ibeg - 1] = val;
                upper = ibeg;
                lower = iend;

                while (upper != lower) {
                        while ((upper < lower) && (val < a[lower-1])) {
                                lower--;
                        }

                        if (upper != lower) {
                                a[upper - 1] = a[lower - 1];
                                b[upper - 1] = b[lower - 1];
                                c[upper - 1] = c[lower - 1];
                                d[upper - 1] = d[lower - 1];
                                indarr[upper - 1] = indarr[lower - 1];
                        }

                        while ((upper < lower) && (val >= a[upper - 1])) {
                                upper++;
                        }

                        if (upper != lower) {
                                a[lower - 1] = a[upper - 1];
                                b[lower - 1] = b[upper - 1];
                                c[lower - 1] = c[upper - 1];
                                d[lower - 1] = d[upper - 1];
                                indarr[lower - 1] = indarr[upper - 1];
                        }
                }

                *midind = upper;

                /* replace TA in position upper and change midind if ta > val */

                if (ta > val)
                        *midind = upper - 1;


                a[upper - 1] = ta;
                b[upper - 1] = tb;
                c[upper - 1] = tc;
                d[upper - 1] = td;
                indarr[upper - 1] = tind;

        } else if (ibeg == iend) {

                if (a[ibeg - 1] < val) 
                        *midind = ibeg;
                else
                        *midind = ibeg - 1;

        } else {

                *midind = ibeg - 1;

        }

        return;

} /* END function cp_partition */



void pc_partition(double *a, double *b, double *c, double *d, double *w, int *indarr,
                  int ibeg, int iend, double val, int *midind)
{
    double ta, tb, tc, td, tw;
    int lower, upper, tind;
    
    
    if (ibeg < iend) {
        /*
         * temporarily stores ibeg entries and set a(ibeg) = val
         * for the partitioning algorithm.
         */
        ta = a[ibeg - 1];
        tb = b[ibeg - 1];
        tc = c[ibeg - 1];
        td = d[ibeg - 1];
        tw = w[ibeg - 1];
        tind = indarr[ibeg - 1];
        a[ibeg - 1] = val;
        upper = ibeg;
        lower = iend;
        
        while (upper != lower) {
            while ((upper < lower) && (val < a[lower-1])) {
                lower--;
            }
            
            if (upper != lower) {
                a[upper - 1] = a[lower - 1];
                b[upper - 1] = b[lower - 1];
                c[upper - 1] = c[lower - 1];
                d[upper - 1] = d[lower - 1];
                w[upper - 1] = w[lower - 1];
                indarr[upper - 1] = indarr[lower - 1];
            }
            
            while ((upper < lower) && (val >= a[upper - 1])) {
                upper++;
            }
            
            if (upper != lower) {
                a[lower - 1] = a[upper - 1];
                b[lower - 1] = b[upper - 1];
                c[lower - 1] = c[upper - 1];
                d[lower - 1] = d[upper - 1];
                w[lower - 1] = w[upper - 1];
                indarr[lower - 1] = indarr[upper - 1];
            }
        }
        
        *midind = upper;
        
        
        /* replace TA in position upper and change midind if ta > val */
        
        if (ta > val)
            *midind = upper - 1;
        
        
        a[upper - 1] = ta;
        b[upper - 1] = tb;
        c[upper - 1] = tc;
        d[upper - 1] = td;
        w[upper - 1] = tw;
        indarr[upper - 1] = tind;
        
    } else if (ibeg == iend) {
        
        if (a[ibeg - 1] < val)
            *midind = ibeg;
        else
            *midind = ibeg - 1;
        
    } else {
        
        *midind = ibeg - 1;
        
    }
    
    return;
    
} /* END function pc_partition */


void pc_partition_8(double *x, double *y, double *z, double *q, double *w, int *orderarr, double xyzmms[6][8],
                    double xl, double yl, double zl, int *numposchild, int max_num_children,
                    double x_mid, double y_mid, double z_mid, int ind[8][2])
{
    int temp_ind;

    *numposchild = 1;
    double critlen = max3(xl, yl, zl) / sqrt(2.0);
    
    int divide_x = 0;
    int divide_y = 0;
    int divide_z = 0;
    
    if (xl >= critlen) divide_x = 1;
    if (yl >= critlen) divide_y = 1;
    if (zl >= critlen) divide_z = 1;
    
    if (max_num_children == 4) {
        if (xl < yl && xl < zl) divide_x = 0;
        if (yl < xl && yl < zl) divide_y = 0;
        if (zl < xl && zl < yl) divide_z = 0;
        
        if (divide_x + divide_y + divide_z == 3) divide_x = 0;
    }
    
    if (max_num_children == 2) {
        if (xl < yl || xl < zl) divide_x = 0;
        if (yl < xl || yl < zl) divide_y = 0;
        if (zl < xl || zl < yl) divide_z = 0;
        
        if (divide_x + divide_y + divide_z == 3) {
            divide_x = 0;
            divide_y = 0;
        }
        
        if (divide_x + divide_y + divide_z == 2) {
            if        (divide_x == 0) {
                divide_y = 0;
            } else if (divide_y == 0) {
                divide_z = 0;
            } else if (divide_z == 0) {
                divide_x = 0;
            }
        }
    }

    if (divide_x) {

        pc_partition(x, y, z, q, w, orderarr, ind[0][0], ind[0][1],
                     x_mid, &temp_ind);

        ind[1][0] = temp_ind + 1;
        ind[1][1] = ind[0][1];
        ind[0][1] = temp_ind;

        for (int i = 0; i < 6; i++)
            xyzmms[i][1] = xyzmms[i][0];

        xyzmms[1][0] = x_mid;
        xyzmms[0][1] = x_mid;
        *numposchild = 2 * *numposchild;

    }

    if (divide_y) {

        for (int i = 0; i < *numposchild; i++) {
            pc_partition(y, x, z, q, w, orderarr, ind[i][0], ind[i][1],
                         y_mid, &temp_ind);

            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (int j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[3][i] = y_mid;
            xyzmms[2][*numposchild + i] = y_mid;
        }

        *numposchild = 2 * *numposchild;

    }

    if (divide_z) {

        for (int i = 0; i < *numposchild; i++) {
            pc_partition(z, x, y, q, w, orderarr, ind[i][0], ind[i][1],
                         z_mid, &temp_ind);

            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (int j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[5][i] = z_mid;
            xyzmms[4][*numposchild + i] = z_mid;
        }

        *numposchild = 2 * *numposchild;

    }

    return;

} /* END of function partition_8 */


void cp_partition_8(double *x, double *y, double *z, double *q, int *orderarr, double xyzmms[6][8],
                    double xl, double yl, double zl, int *numposchild, int max_num_children,
                    double x_mid, double y_mid, double z_mid, int ind[8][2])
{
    int temp_ind;

    *numposchild = 1;
    double critlen = max3(xl, yl, zl) / sqrt(2.0);
    
    int divide_x = 0;
    int divide_y = 0;
    int divide_z = 0;
    
    if (xl >= critlen) divide_x = 1;
    if (yl >= critlen) divide_y = 1;
    if (zl >= critlen) divide_z = 1;
    
    if (max_num_children == 4) {
        if (xl < yl && xl < zl) divide_x = 0;
        if (yl < xl && yl < zl) divide_y = 0;
        if (zl < xl && zl < yl) divide_z = 0;
        
        if (divide_x + divide_y + divide_z == 3) divide_x = 0;
    }
    
    if (max_num_children == 2) {
        if (xl < yl || xl < zl) divide_x = 0;
        if (yl < xl || yl < zl) divide_y = 0;
        if (zl < xl || zl < yl) divide_z = 0;
        
        if (divide_x + divide_y + divide_z == 3) {
            divide_x = 0;
            divide_y = 0;
        }
        
        if (divide_x + divide_y + divide_z == 2) {
            if        (divide_x == 0) {
                divide_y = 0;
            } else if (divide_y == 0) {
                divide_z = 0;
            } else if (divide_z == 0) {
                divide_x = 0;
            }
        }
    }

    if (divide_x) {
        cp_partition(x, y, z, q, orderarr, ind[0][0], ind[0][1],
                     x_mid, &temp_ind);

        ind[1][0] = temp_ind + 1;
        ind[1][1] = ind[0][1];
        ind[0][1] = temp_ind;

        for (int i = 0; i < 6; i++)
            xyzmms[i][1] = xyzmms[i][0];

        xyzmms[1][0] = x_mid;
        xyzmms[0][1] = x_mid;
        *numposchild = 2 * *numposchild;
    }

    if (divide_y) {
        for (int i = 0; i < *numposchild; i++) {
            cp_partition(y, x, z, q, orderarr, ind[i][0], ind[i][1],
                         y_mid, &temp_ind);

            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (int j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[3][i] = y_mid;
            xyzmms[2][*numposchild + i] = y_mid;
        }

        *numposchild = 2 * *numposchild;
    }

    if (divide_z) {
        for (int i = 0; i < *numposchild; i++) {
            cp_partition(z, x, y, q, orderarr, ind[i][0], ind[i][1],
                         z_mid, &temp_ind);

            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (int j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[5][i] = z_mid;
            xyzmms[4][*numposchild + i] = z_mid;
        }

        *numposchild = 2 * *numposchild;

    }

    return;

} /* END of function cp_partition_8 */

