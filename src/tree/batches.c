#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"

#include "../particles/struct_particles.h"

#include "struct_tree.h"
#include "partition.h"
#include "batches.h"


static void Batches_Targets_Fill(struct Tree *batches, int *sizeof_batch_arrays,
                    int maxparnode, double *xyzmm, int *xyzdim, int *xyzind);
                    
static void Batches_Sources_Fill(struct Tree *batches, int *sizeof_batch_arrays, struct Particles *particles,
                    int ibeg, int iend, int maxparnode, double *xyzmm);
                    
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

    make_vector(batches->x_min, init_sizeof_arrays);
    make_vector(batches->y_min, init_sizeof_arrays);
    make_vector(batches->z_min, init_sizeof_arrays);
    
    make_vector(batches->x_max, init_sizeof_arrays);
    make_vector(batches->y_max, init_sizeof_arrays);
    make_vector(batches->z_max, init_sizeof_arrays);

    make_vector(batches->x_dim, init_sizeof_arrays);
    make_vector(batches->y_dim, init_sizeof_arrays);
    make_vector(batches->z_dim, init_sizeof_arrays);

    make_vector(batches->x_low_ind, init_sizeof_arrays);
    make_vector(batches->y_low_ind, init_sizeof_arrays);
    make_vector(batches->z_low_ind, init_sizeof_arrays);

    make_vector(batches->x_high_ind, init_sizeof_arrays);
    make_vector(batches->y_high_ind, init_sizeof_arrays);
    make_vector(batches->z_high_ind, init_sizeof_arrays);
    
    Batches_Sources_Fill(batches, &init_sizeof_arrays, sources, 1, sources->num,
                         run_params->max_per_source_leaf, xyzminmax);

#ifdef OPENACC_ENABLED
    double *source_x = sources->x;
    double *source_y = sources->y;
    double *source_z = sources->z;
    double *source_q = sources->q;
    int source_num = sources->num;
    #pragma acc enter data copyin(source_x[0:source_num], source_y[0:source_num], \
                                  source_z[0:source_num], source_q[0:source_num])
#endif

    return;
}


void Batches_Targets_Construct(struct Tree **batches_addr, struct Particles *targets, struct RunParams *run_params)
{
    double xyzminmax[6];
    int xyzdim[3], xyzind[6];
    
    xyzminmax[0] = minval(targets->x, targets->num);
    xyzminmax[1] = maxval(targets->x, targets->num);
    xyzminmax[2] = minval(targets->y, targets->num);
    xyzminmax[3] = maxval(targets->y, targets->num);
    xyzminmax[4] = minval(targets->z, targets->num);
    xyzminmax[5] = maxval(targets->z, targets->num);
    
    xyzdim[0] = targets->xdim;
    xyzdim[1] = targets->xdim;
    xyzdim[2] = targets->ydim;
    
    xyzind[0] = 0;
    xyzind[1] = targets->xdim-1;
    xyzind[2] = 0;
    xyzind[3] = targets->ydim-1;
    xyzind[4] = 0;
    xyzind[5] = targets->zdim-1;
    
    
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
    
    make_vector(batches->x_min, init_sizeof_arrays);
    make_vector(batches->y_min, init_sizeof_arrays);
    make_vector(batches->z_min, init_sizeof_arrays);
    
    make_vector(batches->x_max, init_sizeof_arrays);
    make_vector(batches->y_max, init_sizeof_arrays);
    make_vector(batches->z_max, init_sizeof_arrays);

    make_vector(batches->x_dim, init_sizeof_arrays);
    make_vector(batches->y_dim, init_sizeof_arrays);
    make_vector(batches->z_dim, init_sizeof_arrays);

    make_vector(batches->x_low_ind, init_sizeof_arrays);
    make_vector(batches->y_low_ind, init_sizeof_arrays);
    make_vector(batches->z_low_ind, init_sizeof_arrays);

    make_vector(batches->x_high_ind, init_sizeof_arrays);
    make_vector(batches->y_high_ind, init_sizeof_arrays);
    make_vector(batches->z_high_ind, init_sizeof_arrays);
    
    Batches_Targets_Fill(batches, &init_sizeof_arrays,
                         run_params->max_per_target_leaf, xyzminmax, xyzdim, xyzind);

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

        make_vector(batches->x_min, length);
        make_vector(batches->y_min, length);
        make_vector(batches->z_min, length);
        
        make_vector(batches->x_max, length);
        make_vector(batches->y_max, length);
        make_vector(batches->z_max, length);

        make_vector(batches->x_dim, length);
        make_vector(batches->y_dim, length);
        make_vector(batches->z_dim, length);

        make_vector(batches->x_low_ind, length);
        make_vector(batches->y_low_ind, length);
        make_vector(batches->z_low_ind, length);

        make_vector(batches->x_high_ind, length);
        make_vector(batches->y_high_ind, length);
        make_vector(batches->z_high_ind, length);
    }

    return;
} /* END of function setup */



void Batches_Free(struct Tree **batches_addr)
{
    struct Tree *batches = *batches_addr;
    
    if (batches != NULL) {
        if (batches->ibeg   != NULL) free_vector(batches->ibeg);
        if (batches->iend   != NULL) free_vector(batches->iend);
        if (batches->numpar != NULL) free_vector(batches->numpar);
        
        if (batches->x_mid  != NULL) free_vector(batches->x_mid);
        if (batches->y_mid  != NULL) free_vector(batches->y_mid);
        if (batches->z_mid  != NULL) free_vector(batches->z_mid);
        if (batches->radius != NULL) free_vector(batches->radius);
        
        if (batches->x_min  != NULL) free_vector(batches->x_min);
        if (batches->y_min  != NULL) free_vector(batches->y_min);
        if (batches->z_min  != NULL) free_vector(batches->z_min);
        
        if (batches->x_max  != NULL) free_vector(batches->x_max);
        if (batches->y_max  != NULL) free_vector(batches->y_max);
        if (batches->z_max  != NULL) free_vector(batches->z_max);
        
        if (batches->x_low_ind != NULL) free_vector(batches->x_low_ind);
        if (batches->y_low_ind != NULL) free_vector(batches->y_low_ind);
        if (batches->z_low_ind != NULL) free_vector(batches->z_low_ind);

        if (batches->x_high_ind != NULL) free_vector(batches->x_high_ind);
        if (batches->y_high_ind != NULL) free_vector(batches->y_high_ind);
        if (batches->z_high_ind != NULL) free_vector(batches->z_high_ind);

        if (batches->x_dim != NULL) free_vector(batches->x_dim);
        if (batches->y_dim != NULL) free_vector(batches->y_dim);
        if (batches->z_dim != NULL) free_vector(batches->z_dim);

        free(batches);
    }
    
    batches = NULL;

    return;
} /* END of function setup */



void Batches_Free_Win(struct Tree **batches_addr)
{
    struct Tree *batches = *batches_addr;
    
    if (batches != NULL) {
        if (batches->ibeg   != NULL) free_vector(batches->ibeg);
        if (batches->iend   != NULL) free_vector(batches->iend);
        if (batches->numpar != NULL) free_vector(batches->numpar);
        
        if (batches->x_mid  != NULL) free_vector(batches->x_mid);
        if (batches->y_mid  != NULL) free_vector(batches->y_mid);
        if (batches->z_mid  != NULL) free_vector(batches->z_mid);
        if (batches->radius != NULL) free_vector(batches->radius);

        if (batches->x_min  != NULL) free_vector(batches->x_min);
        if (batches->y_min  != NULL) free_vector(batches->y_min);
        if (batches->z_min  != NULL) free_vector(batches->z_min);
        
        if (batches->x_max  != NULL) free_vector(batches->x_max);
        if (batches->y_max  != NULL) free_vector(batches->y_max);
        if (batches->z_max  != NULL) free_vector(batches->z_max);
        
        if (batches->x_low_ind != NULL) free_vector(batches->x_low_ind);
        if (batches->y_low_ind != NULL) free_vector(batches->y_low_ind);
        if (batches->z_low_ind != NULL) free_vector(batches->z_low_ind);

        if (batches->x_high_ind != NULL) free_vector(batches->x_high_ind);
        if (batches->y_high_ind != NULL) free_vector(batches->y_high_ind);
        if (batches->z_high_ind != NULL) free_vector(batches->z_high_ind);

        if (batches->x_dim != NULL) free_vector(batches->x_dim);
        if (batches->y_dim != NULL) free_vector(batches->y_dim);
        if (batches->z_dim != NULL) free_vector(batches->z_dim);

        free(batches);
    }
    
    batches = NULL;
    
    return;
} /* END of function setup */



void Batches_Print(struct Tree *batches)
{
    printf("[BaryTree]\n");
    printf("[BaryTree] Batches information: \n");
    printf("[BaryTree]\n");
    printf("[BaryTree]             Cumulative batches: %d\n", batches->numnodes);
    printf("[BaryTree]\n");

    return;
}




//~~~~~~~~~~~~~~~~~~~~~~~~~~//
//~~~LOCAL FUNCTIONS~~~~~~~~//
//~~~~~~~~~~~~~~~~~~~~~~~~~~//

void Batches_Targets_Fill(struct Tree *batches, int *sizeof_batch_arrays,
                          int maxparnode, double *xyzmm, int *xyzdim, int *xyzind)
{
    int numpar = xyzdim[0]*xyzdim[1]*xyzdim[2];
    
    double x_min = xyzmm[0];
    double x_max = xyzmm[1];
    double y_min = xyzmm[2];
    double y_max = xyzmm[3];
    double z_min = xyzmm[4];
    double z_max = xyzmm[5];

    double xl = x_max - x_min;
    double yl = y_max - y_min;
    double zl = z_max - z_min;
    
    double x_mid = (x_max + x_min) / 2.0;
    double y_mid = (y_max + y_min) / 2.0;
    double z_mid = (z_max + z_min) / 2.0;

    double radius = sqrt(xl*xl + yl*yl + zl*zl) / 2.0;

    if (numpar > maxparnode) {
    
        int max_num_children;
    
        if (numpar < 2 * maxparnode) {
            max_num_children = 2;
        } else if (numpar < 4 * maxparnode) {
            max_num_children = 4;
        } else {
            max_num_children = 8;
        }
        
        int xyzdims[3][8], xyzinds[6][8];
        double xyzmms[6][8];
    
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 8; j++)
                xyzdims[i][j] = 0;
    
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 8; j++) {
                xyzmms[i][j] = 0.0;
                xyzinds[i][j] = 0;
            }
        }
    /*
     * IND array holds indices of the eight new subregions.
     */
        xyzmms[0][0] = x_min;
        xyzmms[1][0] = x_max;
        xyzmms[2][0] = y_min;
        xyzmms[3][0] = y_max;
        xyzmms[4][0] = z_min;
        xyzmms[5][0] = z_max;

        xyzdims[0][0] = xyzdim[0];
        xyzdims[1][0] = xyzdim[1];
        xyzdims[2][0] = xyzdim[2];
        
        xyzinds[0][0] = xyzind[0];
        xyzinds[1][0] = xyzind[1];
        xyzinds[2][0] = xyzind[2];
        xyzinds[3][0] = xyzind[3];
        xyzinds[4][0] = xyzind[4];
        xyzinds[5][0] = xyzind[5];

        int numposchild;
        
        cp_partition_8(xyzmms, xyzdims, xyzinds, xl, yl, zl, &numposchild, max_num_children);

        for (int i = 0; i < numposchild; i++) {
            if (xyzinds[0][i] <= xyzinds[1][i] &&
                xyzinds[2][i] <= xyzinds[3][i] &&
                xyzinds[4][i] <= xyzinds[5][i]) {

                double lxyzmm[6];
                int lxyzind[6], lxyzdim[3];

                for (int j = 0; j < 6; j++) {
                    lxyzmm[j] = xyzmms[j][i];
                    lxyzind[j] = xyzinds[j][i];
                }

                for (int j = 0; j < 3; j++) {
                    lxyzdim[j] = xyzdims[j][i];
                }
                               
                Batches_Targets_Fill(batches, sizeof_batch_arrays,
                                     maxparnode, lxyzmm, lxyzdim, lxyzind);
            }
        }

    } else {
        
        if (batches->numnodes >= *sizeof_batch_arrays) {
            (*sizeof_batch_arrays) *= 1.5;
            Batches_ReallocArrays(batches, *sizeof_batch_arrays);
        }
        
        batches->numpar[batches->numnodes] = numpar;
        
        batches->x_mid[batches->numnodes] = x_mid;
        batches->y_mid[batches->numnodes] = y_mid;
        batches->z_mid[batches->numnodes] = z_mid;
        
        batches->radius[batches->numnodes] = radius;
        
        
        batches->x_min[batches->numnodes] = x_min;
        batches->y_min[batches->numnodes] = y_min;
        batches->z_min[batches->numnodes] = z_min;
        
        batches->x_max[batches->numnodes] = x_max;
        batches->y_max[batches->numnodes] = y_max;
        batches->z_max[batches->numnodes] = z_max;
    
        batches->x_dim[batches->numnodes] = xyzdim[0];
        batches->y_dim[batches->numnodes] = xyzdim[1];
        batches->z_dim[batches->numnodes] = xyzdim[2];
    
        batches->x_low_ind[batches->numnodes] = xyzind[0];
        batches->y_low_ind[batches->numnodes] = xyzind[2];
        batches->z_low_ind[batches->numnodes] = xyzind[4];
    
        batches->x_high_ind[batches->numnodes] = xyzind[1];
        batches->y_high_ind[batches->numnodes] = xyzind[3];
        batches->z_high_ind[batches->numnodes] = xyzind[5];
        
        batches->numnodes++;
    }

    return;

} /* end of function create_target_batch */



void Batches_Sources_Fill(struct Tree *batches, int *sizeof_batch_arrays, struct Particles *particles,
                          int ibeg, int iend, int maxparnode, double *xyzmm)
{
    int numposchild;
    
    int ind[8][2];
    double xyzmms[6][8];
    double lxyzmm[6];

    for (int i = 0; i < 8; i++) {
        for (int j = 0; j < 2; j++) {
            ind[i][j] = 0.0;
        }
    }

    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 8; j++) {
            xyzmms[i][j] = 0.0;
        }
    }

    for (int i = 0; i < 6; i++) {
        lxyzmm[i] = 0.0;
    }
    
    int numpar = iend - ibeg + 1;

    double x_min = minval(particles->x + ibeg - 1, numpar);
    double x_max = maxval(particles->x + ibeg - 1, numpar);
    double y_min = minval(particles->y + ibeg - 1, numpar);
    double y_max = maxval(particles->y + ibeg - 1, numpar);
    double z_min = minval(particles->z + ibeg - 1, numpar);
    double z_max = maxval(particles->z + ibeg - 1, numpar);
    
    double xl = x_max - x_min;
    double yl = y_max - y_min;
    double zl = z_max - z_min;
    
    double x_mid = (x_max + x_min) / 2.0;
    double y_mid = (y_max + y_min) / 2.0;
    double z_mid = (z_max + z_min) / 2.0;

    double radius = sqrt(xl*xl + yl*yl + zl*zl) / 2.0;
    
    if (numpar > maxparnode) {
    
        int max_num_children;
    
        if (numpar < 2 * maxparnode) {
            max_num_children = 2;
        } else if (numpar < 4 * maxparnode) {
            max_num_children = 4;
        } else {
            max_num_children = 8;
        }
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

        pc_partition_8(particles->x, particles->y, particles->z, particles->q,
                       particles->order, xyzmms, xl, yl, zl, &numposchild, max_num_children,
                       x_mid, y_mid, z_mid, ind);

        for (int i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {

                for (int j = 0; j < 6; j++)
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

        realloc_vector(batches->x_min, newlength);
        realloc_vector(batches->y_min, newlength);
        realloc_vector(batches->z_min, newlength);

        realloc_vector(batches->x_max, newlength);
        realloc_vector(batches->y_max, newlength);
        realloc_vector(batches->z_max, newlength);

        realloc_vector(batches->x_dim, newlength);
        realloc_vector(batches->y_dim, newlength);
        realloc_vector(batches->z_dim, newlength);

        realloc_vector(batches->x_low_ind, newlength);
        realloc_vector(batches->y_low_ind, newlength);
        realloc_vector(batches->z_low_ind, newlength);

        realloc_vector(batches->x_high_ind, newlength);
        realloc_vector(batches->y_high_ind, newlength);
        realloc_vector(batches->z_high_ind, newlength);
    }
    
    return;
}
