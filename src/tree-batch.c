/*
 *Procedures for Cluster-Particle Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"



void cp_setup_batch(double *x, double *y, double *z,
           int numpars, int batch_size, double *xyzminmax, int **batch_reorder,
           int *batch_num, int ***batch_index, double ***batch_center,
           double **batch_radius)
{
    /* local variables */
    int i;
    int max_batch_num;
    
    /* find bounds of Cartesian box enclosing the particles */
    xyzminmax[0] = minval(x, numpars);
    xyzminmax[1] = maxval(x, numpars);
    xyzminmax[2] = minval(y, numpars);
    xyzminmax[3] = maxval(y, numpars);
    xyzminmax[4] = minval(z, numpars);
    xyzminmax[5] = maxval(z, numpars);
    
    *batch_num = 0;
    max_batch_num = (int)ceil((double)numpars * 8 / batch_size);

    make_vector(*batch_reorder, numpars);
    make_matrix(*batch_index, max_batch_num, 2);
    make_matrix(*batch_center, max_batch_num, 3);
    make_vector(*batch_radius, max_batch_num);

    for (i = 0; i < numpars; i++)
        (*batch_reorder)[i] = i+1;

    return;
    
} /* END of function setup */




void cp_create_batch(struct tnode **p, int ibeg, int iend,
                    double *x, double *y, double *z,
                    int maxparnode, double *xyzmm, int level,
                    int *batch_reorder, int *batch_num,
                    int **batch_index, double **batch_center,
                    double *batch_radius)
{
    /*local variables*/
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
    int i, j, loclev, numposchild, nump;
    
    int ind[8][2];
    double xyzmms[6][8];
    double lxyzmm[6];

    for (i = 0; i < 8; i++) {
        for (j = 0; j < 2; j++) {
            ind[i][j] = 0.0;
        }
    }

    for (i = 0; i < 6; i++) {
        for (j = 0; j < 8; j++) {
            xyzmms[i][j] = 0.0;
        }
    }

    for (i = 0; i < 6; i++) {
        lxyzmm[i] = 0.0;
    }

    (*p) = malloc(sizeof(struct tnode));

    /* set node fields: number of particles, exist_ms, and xyz bounds */
    (*p)->numpar = iend - ibeg + 1;

    (*p)->x_min = minval(x + ibeg - 1, (*p)->numpar);
    (*p)->x_max = maxval(x + ibeg - 1, (*p)->numpar);
    (*p)->y_min = minval(y + ibeg - 1, (*p)->numpar);
    (*p)->y_max = maxval(y + ibeg - 1, (*p)->numpar);
    (*p)->z_min = minval(z + ibeg - 1, (*p)->numpar);
    (*p)->z_max = maxval(z + ibeg - 1, (*p)->numpar);

    /*compute aspect ratio*/
    xl = (*p)->x_max - (*p)->x_min;
    yl = (*p)->y_max - (*p)->y_min;
    zl = (*p)->z_max - (*p)->z_min;
        
    lmax = max3(xl, yl, zl);
    t1 = lmax;
    t2 = min3(xl, yl, zl);


    if (t2 != 0.0)
        (*p)->aspect = t1/t2;
    else
        (*p)->aspect = 0.0;

    
    /*midpoint coordinates, RADIUS and SQRADIUS*/
    (*p)->x_mid = ((*p)->x_max + (*p)->x_min) / 2.0;
    (*p)->y_mid = ((*p)->y_max + (*p)->y_min) / 2.0;
    (*p)->z_mid = ((*p)->z_max + (*p)->z_min) / 2.0;

    t1 = (*p)->x_max - (*p)->x_mid;
    t2 = (*p)->y_max - (*p)->y_mid;
    t3 = (*p)->z_max - (*p)->z_mid;

    (*p)->sqradius = t1*t1 + t2*t2 + t3*t3;
    (*p)->radius = sqrt((*p)->sqradius);

    
    /*set particle limits, tree level of node, and nullify child pointers*/
    (*p)->ibeg = ibeg;
    (*p)->iend = iend;
    (*p)->level = level;

    (*p)->num_children = 0;
    for (i = 0; i < 8; i++)
        (*p)->child[i] = NULL;
    

    if ((*p)->numpar > maxparnode) {
    /*
     * set IND array to 0, and then call PARTITION_8 routine.
     * IND array holds indices of the eight new subregions.
     * Also, setup XYZMMS array in the case that SHRINK = 1.
     */
        xyzmms[0][0] = (*p)->x_min;
        xyzmms[1][0] = (*p)->x_max;
        xyzmms[2][0] = (*p)->y_min;
        xyzmms[3][0] = (*p)->y_max;
        xyzmms[4][0] = (*p)->z_min;
        xyzmms[5][0] = (*p)->z_max;

        ind[0][0] = ibeg;
        ind[0][1] = iend;

        x_mid = (*p)->x_mid;
        y_mid = (*p)->y_mid;
        z_mid = (*p)->z_mid;


        cp_partition_batch(x, y, z, xyzmms, xl, yl, zl,
                           lmax, &numposchild,
                           x_mid, y_mid, z_mid, ind, batch_reorder);

        loclev = level + 1;

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {
                (*p)->num_children = (*p)->num_children + 1;

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];

                cp_create_batch(&((*p)->child[(*p)->num_children - 1]),
                                  ind[i][0], ind[i][1], x, y, z,
                                  maxparnode, lxyzmm, loclev,
                                  batch_reorder, batch_num, batch_index, batch_center,
                                  batch_radius);
            }
        }

    } else {
    
        (*batch_num) += 1;
        
        batch_index[*batch_num-1][0] = (*p)->ibeg;
        batch_index[*batch_num-1][1] = (*p)->iend;
        
        batch_center[*batch_num-1][0] = (*p)->x_mid;
        batch_center[*batch_num-1][1] = (*p)->y_mid;
        batch_center[*batch_num-1][2] = (*p)->z_mid;
        
        batch_radius[*batch_num-1] = (*p)->radius;
    }

    return;

} /* end of function create_tree_n0 */




void cp_partition_batch(double *x, double *y, double *z, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *batch_reorder)
{

    /* local variables */
    int temp_ind, i, j;
    double critlen;

    *numposchild = 1;
    critlen = lmax / sqrt(2.0);

    if (xl >= critlen) {
        cp_partition(x, y, z, batch_reorder, ind[0][0], ind[0][1],
                     x_mid, &temp_ind);

        ind[1][0] = temp_ind + 1;
        ind[1][1] = ind[0][1];
        ind[0][1] = temp_ind;

        for (i = 0; i < 6; i++)
            xyzmms[i][1] = xyzmms[i][0];
        
        xyzmms[1][0] = x_mid;
        xyzmms[0][1] = x_mid;
        *numposchild = 2 * *numposchild;
    }

    if (yl >= critlen) {
        for (i = 0; i < *numposchild; i++) {
            cp_partition(y, x, z, batch_reorder, ind[i][0], ind[i][1],
                         y_mid, &temp_ind);
            
            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[3][i] = y_mid;
            xyzmms[2][*numposchild + i] = y_mid;
        }
        
        *numposchild = 2 * *numposchild;
    }

    if (zl >= critlen) {
        for (i = 0; i < *numposchild; i++) {
            cp_partition(z, x, y, batch_reorder, ind[i][0], ind[i][1],
                         z_mid, &temp_ind);
            
            ind[*numposchild + i][0] = temp_ind + 1;
            ind[*numposchild + i][1] = ind[i][1];
            ind[i][1] = temp_ind;

            for (j = 0; j < 6; j++)
                xyzmms[j][*numposchild + i] = xyzmms[j][i];

            xyzmms[5][i] = z_mid;
            xyzmms[4][*numposchild + i] = z_mid;
        }
        
        *numposchild = 2 * *numposchild;

    }

    return;

} /* END of function cp_partition_8 */

