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
#include "particles.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"



void cp_partition_batch(double *x, double *y, double *z, double *q, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *batch_reorder);

void pc_partition_batch(double *x, double *y, double *z, double *q, double *w, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2],
                    int *batch_reorder);



void setup_batch(struct batch **batches, double *batch_lim,
                 struct particles *particles, int batch_size)
{
    /* local variables */
    int i;
    int max_batch_num;
    
    /* find bounds of Cartesian box enclosing the particles */
    batch_lim[0] = minval(particles->x, particles->num);
    batch_lim[1] = maxval(particles->x, particles->num);
    batch_lim[2] = minval(particles->y, particles->num);
    batch_lim[3] = maxval(particles->y, particles->num);
    batch_lim[4] = minval(particles->z, particles->num);
    batch_lim[5] = maxval(particles->z, particles->num);
    
    (*batches) = malloc(sizeof(struct batch));
    (*batches)->numnodes = 0;
    
    max_batch_num = 2*(int)ceil((double)particles->num * 8 / batch_size); // this needs improvement.  I set to 2* to stop it from crashing.

    make_vector((*batches)->reorder, particles->num);

    //make_matrix((*batches)->index, max_batch_num, 4);
    make_vector((*batches)->ibeg, max_batch_num);
    make_vector((*batches)->iend, max_batch_num);
    make_vector((*batches)->numpar, max_batch_num);
    make_vector((*batches)->numApprox, max_batch_num);
    make_vector((*batches)->numDirect, max_batch_num);

    //make_matrix((*batches)->center, max_batch_num, 3);
    make_vector((*batches)->x_mid, max_batch_num);
    make_vector((*batches)->y_mid, max_batch_num);
    make_vector((*batches)->z_mid, max_batch_num);

    make_vector((*batches)->radius, max_batch_num);

    for (i = 0; i < particles->num; i++)
        (*batches)->reorder[i] = i+1;

    return;
    
} /* END of function setup */




void create_target_batch(struct batch *batches, struct particles *particles,
                         int ibeg, int iend, int maxparnode, double *xyzmm)
{
	int rank; int numProcs;	int ierr;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
//	printf("Entered create_target_batch.\n");
//	fflush(stdout);
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
//    if (rank==0){
//    	printf("numpar = %i\n", numpar);
//    }
    x_min = minval(particles->x + ibeg - 1, numpar);
    x_max = maxval(particles->x + ibeg - 1, numpar);
    y_min = minval(particles->y + ibeg - 1, numpar);
    y_max = maxval(particles->y + ibeg - 1, numpar);
    z_min = minval(particles->z + ibeg - 1, numpar);
    z_max = maxval(particles->z + ibeg - 1, numpar);
    
//    if (rank==1){
//		printf("x_min, x_max = %f, %f\n", x_min, x_max);
//	}

//    printf("Got here.\n");
    /*compute aspect ratio*/
    xl = x_max - x_min;
    yl = y_max - y_min;
    zl = z_max - z_min;
    
    lmax = max3(xl, yl, zl);
    
    /*midpoint coordinates, RADIUS and SQRADIUS*/
    x_mid = (x_max + x_min) / 2.0;
    y_mid = (y_max + y_min) / 2.0;
    z_mid = (z_max + z_min) / 2.0;



//    printf("xmid, ymid, zmid = %1.2f, %1.2f, %1.2f.\n", x_mid,y_mid,z_mid);
//	printf("xmin, xmax, ymin, ymax, zmin, zmax = %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f.\n", x_min,x_max,y_min,y_max,z_min,z_max);


    t1 = x_max - x_mid;
    t2 = y_max - y_mid;
    t3 = z_max - z_mid;

    sqradius = t1*t1 + t2*t2 + t3*t3;
    radius = sqrt(sqradius);
    
//    printf("Set bounds and other variables.\n");

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



        cp_partition_batch(particles->x, particles->y, particles->z, particles->q,
                           xyzmms, xl, yl, zl, lmax, &numposchild,
                           x_mid, y_mid, z_mid, ind, batches->reorder);

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];
                
                create_target_batch(batches, particles, ind[i][0], ind[i][1],
                                    maxparnode, lxyzmm);

//                printf("Creating child %i.\n", i);
            }
        }

    } else {

        batches->numnodes += 1;
//    	printf("Rank %i Increasing batch number to %i.\n",rank, batches->num);
        
//    	if (batches->num==1){
//    		printf("ibeg %i\n", ibeg);
//    		printf("iend %i\n", iend);
//    		printf("x_mid %f\n", x_mid);
//    		printf("y_mid %f\n", y_mid);
//    		printf("z_mid %f\n", z_mid);
//    		printf("radius %f\n", radius);
//
//    	}
        batches->ibeg[batches->numnodes-1] = ibeg;
        batches->iend[batches->numnodes-1] = iend;
        batches->numpar[batches->numnodes-1] = iend-ibeg+1;
        batches->numApprox[batches->numnodes-1] = 0;
        batches->numDirect[batches->numnodes-1] = 0;
        
        batches->x_mid[batches->numnodes-1] = x_mid;
        batches->y_mid[batches->numnodes-1] = y_mid;
        batches->z_mid[batches->numnodes-1] = z_mid;
        
        batches->radius[batches->numnodes-1] = radius;
//    	printf("Finished filling in batch details.\n");
    }

    return;

} /* end of function create_target_batch */




void create_source_batch(struct batch *batches, struct particles *particles,
                         int ibeg, int iend, int maxparnode, double *xyzmm)
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

        pc_partition_batch(particles->x, particles->y, particles->z, particles->q, particles->w,
                           xyzmms, xl, yl, zl, lmax, &numposchild,
                           x_mid, y_mid, z_mid, ind, batches->reorder);

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];
                
                create_source_batch(batches, particles, ind[i][0], ind[i][1],
                                    maxparnode, lxyzmm);
            }
        }

    } else {
    
        batches->numnodes += 1;

        batches->ibeg[batches->numnodes-1] = ibeg;
        batches->iend[batches->numnodes-1] = iend;
        batches->numpar[batches->numnodes-1] = iend-ibeg+1;
        batches->numApprox[batches->numnodes-1] = 0;
        batches->numDirect[batches->numnodes-1] = 0;
        
        batches->x_mid[batches->numnodes-1] = x_mid;
        batches->y_mid[batches->numnodes-1] = y_mid;
        batches->z_mid[batches->numnodes-1] = z_mid;
        
        batches->radius[batches->numnodes-1] = radius;
    }

    return;

} /* end of function create_source_batch */




void cp_partition_batch(double *x, double *y, double *z, double *q, double xyzmms[6][8],
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
        cp_partition(x, y, z, q, batch_reorder, ind[0][0], ind[0][1],
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
            cp_partition(y, x, z, q, batch_reorder, ind[i][0], ind[i][1],
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
            cp_partition(z, x, y, q, batch_reorder, ind[i][0], ind[i][1],
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




void pc_partition_batch(double *x, double *y, double *z, double *q, double *w, double xyzmms[6][8],
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
        pc_partition(x, y, z, q, w, batch_reorder, ind[0][0], ind[0][1],
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
            pc_partition(y, x, z, q, w, batch_reorder, ind[i][0], ind[i][1],
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
            pc_partition(z, x, y, q, w, batch_reorder, ind[i][0], ind[i][1],
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




void reorder_targets_and_potential(struct particles *targets, double *tEn,
                     int *reorder, int numpars)
{
    double *temp_energy;
    make_vector(temp_energy, numpars);
    for (int i = 0; i < numpars; i++) temp_energy[i] = tEn[i];
    for (int i = 0; i < numpars; i++) tEn[reorder[i]-1] = temp_energy[i];
    free_vector(temp_energy);

    double *temp_x;
    make_vector(temp_x, numpars);
    for (int i = 0; i < numpars; i++) temp_x[i] = targets->x[i];
    for (int i = 0; i < numpars; i++) targets->x[reorder[i]-1] = temp_x[i];
    free_vector(temp_x);

    double *temp_y;
    make_vector(temp_y, numpars);
    for (int i = 0; i < numpars; i++) temp_y[i] = targets->y[i];
    for (int i = 0; i < numpars; i++) targets->y[reorder[i]-1] = temp_y[i];
    free_vector(temp_y);

    double *temp_z;
    make_vector(temp_z, numpars);
    for (int i = 0; i < numpars; i++) temp_z[i] = targets->z[i];
    for (int i = 0; i < numpars; i++) targets->z[reorder[i]-1] = temp_z[i];
    free_vector(temp_z);

    double *temp_q;
    make_vector(temp_q, numpars);
    for (int i = 0; i < numpars; i++) temp_q[i] = targets->q[i];
    for (int i = 0; i < numpars; i++) targets->q[reorder[i]-1] = temp_q[i];
    free_vector(temp_q);

}
