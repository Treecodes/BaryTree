/*
 *Procedures for Cluster-Particle Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "array.h"
#include "globvars.h"
#include "struct_nodes.h"
#include "struct_particles.h"
#include "tools.h"

#include "partition.h"

#include "batches.h"



static void cp_partition_batch(double *x, double *y, double *z, double *q, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *reorder);

static void pc_partition_batch(double *x, double *y, double *z, double *q, double *w, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *reorder);



void Batches_Alloc(struct tnode_array **new_batches, double *batch_lim,
                   struct particles *particles, int batch_size)
{
    
    /* find bounds of Cartesian box enclosing the particles */
    batch_lim[0] = minval(particles->x, particles->num);
    batch_lim[1] = maxval(particles->x, particles->num);
    batch_lim[2] = minval(particles->y, particles->num);
    batch_lim[3] = maxval(particles->y, particles->num);
    batch_lim[4] = minval(particles->z, particles->num);
    batch_lim[5] = maxval(particles->z, particles->num);
    
    *new_batches = malloc(sizeof(struct tnode_array));
    struct tnode_array *batches = *new_batches;

    batches->numnodes = 0;
    // this still needs improvement!
    int max_batch_num = 2*(int)ceil((double)particles->num * 8 / batch_size);

    make_vector(batches->ibeg, max_batch_num);
    make_vector(batches->iend, max_batch_num);
    make_vector(batches->numpar, max_batch_num);
    make_vector(batches->numApprox, max_batch_num);
    make_vector(batches->numDirect, max_batch_num);

    make_vector(batches->x_mid, max_batch_num);
    make_vector(batches->y_mid, max_batch_num);
    make_vector(batches->z_mid, max_batch_num);
    make_vector(batches->radius, max_batch_num);

    return;
    
} /* END of function setup */


void Batches_Free(struct tnode_array *batches)
{

    if (batches != NULL) {
        free_vector(batches->iend);
        free_vector(batches->ibeg);
        free_vector(batches->numpar);
        free_vector(batches->numApprox);
        free_vector(batches->numDirect);
        free_vector(batches->x_mid);
        free_vector(batches->y_mid);
        free_vector(batches->z_mid);
        free_vector(batches->radius);
        free(batches);
    }

    return;

} /* END of function setup */




void Batches_CreateTargetBatches(struct tnode_array *batches, struct particles *particles,
                                 int ibeg, int iend, int maxparnode, double *xyzmm)
{
	int rank; int numProcs;	int ierr;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

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
    
    numpar = iend - ibeg + 1;

    x_min = minval(particles->x + ibeg - 1, numpar);
    x_max = maxval(particles->x + ibeg - 1, numpar);
    y_min = minval(particles->y + ibeg - 1, numpar);
    y_max = maxval(particles->y + ibeg - 1, numpar);
    z_min = minval(particles->z + ibeg - 1, numpar);
    z_max = maxval(particles->z + ibeg - 1, numpar);

    /*compute aspect ratio*/
    xl = x_max - x_min;
    yl = y_max - y_min;
    zl = z_max - z_min;
    
    lmax = max3(xl, yl, zl);
    
    x_mid = (x_max + x_min) / 2.0;
    y_mid = (y_max + y_min) / 2.0;
    z_mid = (z_max + z_min) / 2.0;

    t1 = x_max - x_mid;
    t2 = y_max - y_mid;
    t3 = z_max - z_mid;

    sqradius = t1*t1 + t2*t2 + t3*t3;
    radius = sqrt(sqradius);
    

    if (numpar > maxparnode) {
    /*
     * IND array holds indices of the eight new subregions.
     */
        xyzmms[0][0] = x_min;
        xyzmms[1][0] = x_max;
        xyzmms[2][0] = y_min;
        xyzmms[3][0] = y_max;
        xyzmms[4][0] = z_min;
        xyzmms[5][0] = z_max;

        ind[0][0] = ibeg;
        ind[0][1] = iend;

        cp_partition_batch(particles->x, particles->y, particles->z, particles->q,
                           xyzmms, xl, yl, zl, lmax, &numposchild,
                           x_mid, y_mid, z_mid, ind, particles->order);

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];
                
                Batches_CreateTargetBatches(batches, particles, ind[i][0], ind[i][1],
                                    maxparnode, lxyzmm);

            }
        }

    } else {

        batches->numnodes += 1;
        
        batches->ibeg[batches->numnodes-1] = ibeg;
        batches->iend[batches->numnodes-1] = iend;
        
        batches->x_mid[batches->numnodes-1] = x_mid;
        batches->y_mid[batches->numnodes-1] = y_mid;
        batches->z_mid[batches->numnodes-1] = z_mid;
        
        batches->radius[batches->numnodes-1] = radius;
    }

    return;

} /* end of function create_target_batch */




void Batches_CreateSourceBatches(struct tnode_array *batches, struct particles *particles,
                                 int ibeg, int iend, int maxparnode, double *xyzmm)
{
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
    
    numpar = iend - ibeg + 1;

    x_min = minval(particles->x + ibeg - 1, numpar);
    x_max = maxval(particles->x + ibeg - 1, numpar);
    y_min = minval(particles->y + ibeg - 1, numpar);
    y_max = maxval(particles->y + ibeg - 1, numpar);
    z_min = minval(particles->z + ibeg - 1, numpar);
    z_max = maxval(particles->z + ibeg - 1, numpar);
    
    /*compute aspect ratio*/
    xl = x_max - x_min;
    yl = y_max - y_min;
    zl = z_max - z_min;
    
    lmax = max3(xl, yl, zl);
    
    x_mid = (x_max + x_min) / 2.0;
    y_mid = (y_max + y_min) / 2.0;
    z_mid = (z_max + z_min) / 2.0;

    t1 = x_max - x_mid;
    t2 = y_max - y_mid;
    t3 = z_max - z_mid;

    sqradius = t1*t1 + t2*t2 + t3*t3;
    radius = sqrt(sqradius);
    
    if (numpar > maxparnode) {
    /*
     * IND array holds indices of the eight new subregions.
     */
        xyzmms[0][0] = x_min;
        xyzmms[1][0] = x_max;
        xyzmms[2][0] = y_min;
        xyzmms[3][0] = y_max;
        xyzmms[4][0] = z_min;
        xyzmms[5][0] = z_max;

        ind[0][0] = ibeg;
        ind[0][1] = iend;

        pc_partition_batch(particles->x, particles->y, particles->z, particles->q, particles->w,
                           xyzmms, xl, yl, zl, lmax, &numposchild,
                           x_mid, y_mid, z_mid, ind, particles->order);

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];
                
                Batches_CreateSourceBatches(batches, particles, ind[i][0], ind[i][1],
                                          maxparnode, lxyzmm);
            }
        }

    } else {
    
        batches->numnodes += 1;

        batches->ibeg[batches->numnodes-1] = ibeg;
        batches->iend[batches->numnodes-1] = iend;
        
        batches->x_mid[batches->numnodes-1] = x_mid;
        batches->y_mid[batches->numnodes-1] = y_mid;
        batches->z_mid[batches->numnodes-1] = z_mid;
        
        batches->radius[batches->numnodes-1] = radius;
    }

    return;

} /* end of function create_source_batch */



/*********************************/
/******* LOCAL FUNCTIONS *********/
/*********************************/


static void cp_partition_batch(double *x, double *y, double *z, double *q, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *reorder)
{
    int temp_ind, i, j;
    double critlen;

    *numposchild = 1;
    critlen = lmax / sqrt(2.0);

    if (xl >= critlen) {
        cp_partition(x, y, z, q, reorder, ind[0][0], ind[0][1],
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
            cp_partition(y, x, z, q, reorder, ind[i][0], ind[i][1],
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
            cp_partition(z, x, y, q, reorder, ind[i][0], ind[i][1],
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

} /* END of function cp_partition_batch */




static void pc_partition_batch(double *x, double *y, double *z, double *q, double *w, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *reorder)
{
    int temp_ind, i, j;
    double critlen;

    *numposchild = 1;
    critlen = lmax / sqrt(2.0);

    if (xl >= critlen) {
        pc_partition(x, y, z, q, w, reorder, ind[0][0], ind[0][1],
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
            pc_partition(y, x, z, q, w, reorder, ind[i][0], ind[i][1],
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
            pc_partition(z, x, y, q, w, reorder, ind[i][0], ind[i][1],
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

} /* END of function pc_partition_batch */
