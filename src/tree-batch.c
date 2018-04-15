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
#include "batch.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"



void setup_batch(struct batch **batches, double *batch_lim,
                 double *x, double *y, double *z, int numpars, int batch_size)
{
    /* local variables */
    int i;
    int max_batch_num;
    
    /* find bounds of Cartesian box enclosing the particles */
    batch_lim[0] = minval(x, numpars);
    batch_lim[1] = maxval(x, numpars);
    batch_lim[2] = minval(y, numpars);
    batch_lim[3] = maxval(y, numpars);
    batch_lim[4] = minval(z, numpars);
    batch_lim[5] = maxval(z, numpars);
    
    (*batches) = malloc(sizeof(struct batch));
    (*batches)->num = 0;
    
    max_batch_num = (int)ceil((double)numpars * 8 / batch_size);

    make_vector((*batches)->reorder, numpars);
    make_matrix((*batches)->index, max_batch_num, 2);
    make_matrix((*batches)->center, max_batch_num, 3);
    make_vector((*batches)->radius, max_batch_num);

    for (i = 0; i < numpars; i++)
        (*batches)->reorder[i] = i+1;

    return;
    
} /* END of function setup */




void cp_create_batch(struct batch *batches,
                     int ibeg, int iend, double *x, double *y, double *z,
                     int maxparnode, double *xyzmm)
{
    /*local variables*/
    double x_min, x_max, y_min, y_max, z_min, z_max;
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
    double sqradius, radius;
    int i, j, numposchild, numpar;
    
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
    
    /* set node fields: number of particles, exist_ms, and xyz bounds */
    numpar = iend - ibeg + 1;

    x_min = minval(x + ibeg - 1, numpar);
    x_max = maxval(x + ibeg - 1, numpar);
    y_min = minval(y + ibeg - 1, numpar);
    y_max = maxval(y + ibeg - 1, numpar);
    z_min = minval(z + ibeg - 1, numpar);
    z_max = maxval(z + ibeg - 1, numpar);
    
    /*compute aspect ratio*/
    xl = x_max - x_min;
    yl = y_max - y_min;
    zl = z_max - z_min;
    
    lmax = max3(xl, yl, zl);
    
    /*midpoint coordinates, RADIUS and SQRADIUS*/
    x_mid = (x_max + x_min) / 2.0;
    y_mid = (y_max + y_min) / 2.0;
    z_mid = (z_max + z_min) / 2.0;

    t1 = x_max - x_mid;
    t2 = y_max - y_mid;
    t3 = z_max - z_mid;

    sqradius = t1*t1 + t2*t2 + t3*t3;
    radius = sqrt(sqradius);

    
    /*set particle limits, tree level of node, and nullify child pointers*/

    if (numpar > maxparnode) {
    /*
     * set IND array to 0, and then call PARTITION_8 routine.
     * IND array holds indices of the eight new subregions.
     * Also, setup XYZMMS array in the case that SHRINK = 1.
     */
        xyzmms[0][0] = x_min;
        xyzmms[1][0] = x_max;
        xyzmms[2][0] = y_min;
        xyzmms[3][0] = y_max;
        xyzmms[4][0] = z_min;
        xyzmms[5][0] = z_max;

        ind[0][0] = ibeg;
        ind[0][1] = iend;

        cp_partition_batch(x, y, z, xyzmms, xl, yl, zl,
                           lmax, &numposchild,
                           x_mid, y_mid, z_mid, ind, batches->reorder);

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];

                cp_create_batch(batches, ind[i][0], ind[i][1], x, y, z,
                                maxparnode, lxyzmm);
            }
        }

    } else {
    
        batches->num += 1;
        
        batches->index[batches->num-1][0] = ibeg;
        batches->index[batches->num-1][1] = iend;
        
        batches->center[batches->num-1][0] = x_mid;
        batches->center[batches->num-1][1] = y_mid;
        batches->center[batches->num-1][2] = z_mid;
        
        batches->radius[batches->num-1] = radius;
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


void reorder_energies(int *batch_reorder, int numpars, double *tEn)
{
    int i;
    double *temp_energy;
    
    make_vector(temp_energy, numpars);
    
    for (i = 0; i < numpars; i++) {
        temp_energy[i] = tEn[i];
    }
    
    for (i = 0; i < numpars; i++) {
        tEn[batch_reorder[i]-1] = temp_energy[i];
    }
    
    free_vector(temp_energy);

}