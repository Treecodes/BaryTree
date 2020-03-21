#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <mpi.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"

#include "../particles/struct_particles.h"

#include "struct_tree.h"
#include "partition.h"
#include "batches.h"


static void Batches_Targets_Fill(struct Tree *batches, int *sizeof_batch_arrays, struct Particles *particles,
                    int ibeg, int iend, int maxparnode, double *xyzmm);

static void cp_partition_batch(double *x, double *y, double *z, double *q, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2], int *reorder);
                    
static void Batches_Sources_Fill(struct Tree *batches, int *sizeof_batch_arrays, struct Particles *particles,
                    int ibeg, int iend, int maxparnode, double *xyzmm);
                    
static void pc_partition_batch(double *x, double *y, double *z, double *q, double *w, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2], int *reorder);
                    
static void Batches_ReallocArrays(struct Tree *batches, int newlength);
                    
                    
void Batches_Sources_Construct(struct Tree **batches_addr, struct Particles *sources, struct RunParams *run_params)
{
    double xyzminmax[6];
    
    xyzminmax[0] = minval(sources->x, sources->num);
    xyzminmax[1] = maxval(sources->x, sources->num);
    xyzminmax[2] = minval(sources->y, sources->num);
    xyzminmax[3] = maxval(sources->y, sources->num);
    xyzminmax[4] = minval(sources->z, sources->num);
    xyzminmax[5] = maxval(sources->z, sources->num);
    
    *batches_addr = malloc(sizeof(struct Tree));
    struct Tree *batches = *batches_addr;
    
    batches->numnodes = 0;

    int init_sizeof_arrays = 50;
    
    make_vector(batches->ibeg, init_sizeof_arrays);
    make_vector(batches->iend, init_sizeof_arrays);
    make_vector(batches->numpar, init_sizeof_arrays);

    make_vector(batches->x_mid, init_sizeof_arrays);
    make_vector(batches->y_mid, init_sizeof_arrays);
    make_vector(batches->z_mid, init_sizeof_arrays);
    make_vector(batches->radius, init_sizeof_arrays);
    
    Batches_Sources_Fill(batches, &init_sizeof_arrays, sources, 1, sources->num,
                         run_params->max_per_source_leaf, xyzminmax);

    return;
}


void Batches_Targets_Construct(struct Tree **batches_addr, struct Particles *targets, struct RunParams *run_params)
{
    double xyzminmax[6];
    
    xyzminmax[0] = minval(targets->x, targets->num);
    xyzminmax[1] = maxval(targets->x, targets->num);
    xyzminmax[2] = minval(targets->y, targets->num);
    xyzminmax[3] = maxval(targets->y, targets->num);
    xyzminmax[4] = minval(targets->z, targets->num);
    xyzminmax[5] = maxval(targets->z, targets->num);
    
    *batches_addr = malloc(sizeof(struct Tree));
    struct Tree *batches = *batches_addr;
    
    batches->numnodes = 0;

    int init_sizeof_arrays = 50;
    
    make_vector(batches->ibeg, init_sizeof_arrays);
    make_vector(batches->iend, init_sizeof_arrays);
    make_vector(batches->numpar, init_sizeof_arrays);

    make_vector(batches->x_mid, init_sizeof_arrays);
    make_vector(batches->y_mid, init_sizeof_arrays);
    make_vector(batches->z_mid, init_sizeof_arrays);
    make_vector(batches->radius, init_sizeof_arrays);
    
    Batches_Targets_Fill(batches, &init_sizeof_arrays, targets, 1, targets->num,
                         run_params->max_per_target_leaf, xyzminmax);

    return;
}



void Batches_Alloc(struct Tree **batches_addr, int length)
{
    *batches_addr = malloc(sizeof(struct Tree));
    struct Tree *batches = *batches_addr;

    batches->numnodes = length;
    
    if (length > 0) {
        make_vector(batches->ibeg, length);
        make_vector(batches->iend, length);
        make_vector(batches->numpar, length);

        make_vector(batches->x_mid, length);
        make_vector(batches->y_mid, length);
        make_vector(batches->z_mid, length);
        make_vector(batches->radius, length);
    }

    return;
} /* END of function setup */



void Batches_Free(struct Tree *batches)
{
    if (batches != NULL) {
        if (batches->ibeg   != NULL) free_vector(batches->ibeg);
        if (batches->iend   != NULL) free_vector(batches->iend);
        if (batches->numpar != NULL) free_vector(batches->numpar);
        
        if (batches->x_mid  != NULL) free_vector(batches->x_mid);
        if (batches->y_mid  != NULL) free_vector(batches->y_mid);
        if (batches->z_mid  != NULL) free_vector(batches->z_mid);
        if (batches->radius != NULL) free_vector(batches->radius);
        free(batches);
    }

    return;
} /* END of function setup */



void Batches_Free_Win(struct Tree *batches)
{
    if (batches != NULL) {
        if (batches->ibeg   != NULL) MPI_Free_mem(batches->ibeg);
        if (batches->iend   != NULL) MPI_Free_mem(batches->iend);
        if (batches->numpar != NULL) MPI_Free_mem(batches->numpar);
        
        if (batches->x_mid  != NULL) MPI_Free_mem(batches->x_mid);
        if (batches->y_mid  != NULL) MPI_Free_mem(batches->y_mid);
        if (batches->z_mid  != NULL) MPI_Free_mem(batches->z_mid);
        if (batches->radius != NULL) MPI_Free_mem(batches->radius);
        free(batches);
    }

    return;
} /* END of function setup */



/*********************************/
/******* LOCAL FUNCTIONS *********/
/*********************************/

void Batches_Targets_Fill(struct Tree *batches, int *sizeof_batch_arrays, struct Particles *particles,
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

        cp_partition_batch(particles->x, particles->y, particles->z, particles->q,
                           xyzmms, xl, yl, zl, lmax, &numposchild,
                           x_mid, y_mid, z_mid, ind, particles->order);

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];
                
                Batches_Targets_Fill(batches, sizeof_batch_arrays, particles, ind[i][0], ind[i][1],
                                    maxparnode, lxyzmm);

            }
        }

    } else {
        
        if (batches->numnodes >= *sizeof_batch_arrays) {
            (*sizeof_batch_arrays) *= 1.5;
            Batches_ReallocArrays(batches, *sizeof_batch_arrays);
        }
        
        batches->ibeg[batches->numnodes] = ibeg;
        batches->iend[batches->numnodes] = iend;
        batches->numpar[batches->numnodes] = iend-ibeg+1;
        
        batches->x_mid[batches->numnodes] = x_mid;
        batches->y_mid[batches->numnodes] = y_mid;
        batches->z_mid[batches->numnodes] = z_mid;
        
        batches->radius[batches->numnodes] = radius;
        
        batches->numnodes++;
    }

    return;

} /* end of function create_target_batch */



void Batches_Sources_Fill(struct Tree *batches, int *sizeof_batch_arrays, struct Particles *particles,
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
                
                Batches_Sources_Fill(batches, sizeof_batch_arrays, particles, ind[i][0], ind[i][1],
                                     maxparnode, lxyzmm);
            }
        }

    } else {
    
        if (batches->numnodes >= *sizeof_batch_arrays) {
            (*sizeof_batch_arrays) *= 1.5;
            Batches_ReallocArrays(batches, *sizeof_batch_arrays);
        }
        
        batches->ibeg[batches->numnodes] = ibeg;
        batches->iend[batches->numnodes] = iend;
        batches->numpar[batches->numnodes] = iend-ibeg+1;
        
        batches->x_mid[batches->numnodes] = x_mid;
        batches->y_mid[batches->numnodes] = y_mid;
        batches->z_mid[batches->numnodes] = z_mid;
        
        batches->radius[batches->numnodes] = radius;
        
        batches->numnodes++;
    }

    return;

} /* end of function create_source_batch */



static void Batches_ReallocArrays(struct Tree *batches, int newlength)
{
    if (newlength > 0) {
        realloc_vector(batches->ibeg, newlength);
        realloc_vector(batches->iend, newlength);
        realloc_vector(batches->numpar, newlength);

        realloc_vector(batches->x_mid, newlength);
        realloc_vector(batches->y_mid, newlength);
        realloc_vector(batches->z_mid, newlength);
        realloc_vector(batches->radius, newlength);
    }
    
    return;
}



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
