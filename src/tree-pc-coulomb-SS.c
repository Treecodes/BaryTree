/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"
#include "mkl.h"



void pc_treecode_coulomb_SS(struct tnode *p, struct batch *batches,
                     struct particles *sources, struct particles *targets,
                     double kappaSq, double *tpeng, double *EnP)
{
    /* local variables */
    int i, j;
    
    for (i = 0; i < targets->num; i++)
        EnP[i] = 0.0;
    
    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < p->num_children; j++) {
            compute_pc_coulomb_SS(p->child[j],
                batches->index[i], batches->center[i], batches->radius[i],
                sources->x, sources->y, sources->z, sources->q, sources->w,
                targets->x, targets->y, targets->z, targets->q, kappaSq, EnP);
        }
    }
    
    *tpeng = sum(EnP, targets->num);
    
    return;
    
} /* END of function pc_treecode */




void compute_pc_coulomb_SS(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double kappaSq, double *EnP)
{
//	printf("Entering compute_cp.  Batch start: %d.  Batch end: %d.\n", batch_ind[0]-1, batch_ind[1]);
    /* local variables */
    double dist;
    double tx, ty, tz;
    int i, j;
    double *temp_i, *temp_j, *temp_k;

    /* determine DIST for MAC test */
    tx = batch_mid[0] - p->x_mid;
    ty = batch_mid[1] - p->y_mid;
    tz = batch_mid[2] - p->z_mid;
    dist = sqrt(tx*tx + ty*ty + tz*tz);


    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00)) {
    /*
     * If MAC is accepted and there is more than n0 particles
     * in the box, use the expansion for the approximation.
     */
//    	printf("MAC accepted, performing particle-cluster approximation...\n");

        make_vector(temp_i, torderlim);
        make_vector(temp_j, torderlim);
        make_vector(temp_k, torderlim);

//        printf("Does p->exist_ms? %d\n", p->exist_ms);
        if (p->exist_ms == 0) {
            make_vector(p->ms, (torderlim)*(torderlim)*(torderlim));
            make_vector(p->ms2, (torderlim)*(torderlim)*(torderlim));
            make_vector(p->tx, torderlim);
            make_vector(p->ty, torderlim);
            make_vector(p->tz, torderlim);
            
//            printf("Allocated vectors for ms, tx, ty, tz.\n");


            for (i = 0; i < (torderlim)*(torderlim)*(torderlim); i++){
                p->ms[i] = 0.0;
                p->ms2[i] = 0.0;
            }
//            printf("Zeroed out p->ms \n");

            pc_comp_ms(p, xS, yS, zS, qS, wS);
//            pc_comp_weights(p);
            p->exist_ms = 1;
//            printf("Created 'moments' for node.\n");
        }
        





        // Allocate the matrices and vectors
//        printf("Working on batch from %d to %d\n\n", batch_ind[0] - 1,batch_ind[1]);
        int numberOfTargets = batch_ind[1] - batch_ind[0] + 1;
        int numberOfInterpolationPoints = torderlim*torderlim*torderlim;

        double *kernelMatrix 		= (double *)mkl_malloc(numberOfTargets * numberOfInterpolationPoints * sizeof(double),64);
//        double *kernelMatrix 		= (double *)malloc(numberOfTargets * numberOfInterpolationPoints * sizeof(double));
        double *interactionResult 	= (double *)mkl_malloc(numberOfTargets * sizeof(double),64);
//        double *interactionResult 	= (double *)malloc(numberOfTargets * sizeof(double));


        double *interpolationX = (double *)malloc(numberOfInterpolationPoints * sizeof(double));
        double *interpolationY = (double *)malloc(numberOfInterpolationPoints * sizeof(double));
        double *interpolationZ = (double *)malloc(numberOfInterpolationPoints * sizeof(double));
        double *Weights1	   = (double *)mkl_malloc(numberOfInterpolationPoints * sizeof(double),64);
        double *Weights2	   = (double *)mkl_malloc(numberOfInterpolationPoints * sizeof(double),64);



        // Fill in the interpolation point coordinate vectors.  Not necessary, but helps modularize the next step, filling the kernel matrix.

        int kk = -1;
        int k1, k2, k3;
		for (k3 = 0; k3 < torderlim; k3++) {
			for (k2 = 0; k2 < torderlim; k2++) {
				for (k1 = 0; k1 < torderlim; k1++) {
					kk++;
					interpolationX[kk] = p->tx[k1];
					interpolationY[kk] = p->ty[k2];
					interpolationZ[kk] = p->tz[k3];
					Weights1[kk] = p->ms[kk];
					Weights2[kk] = p->ms2[kk];

				}
			}
		}
//		printf("Filled in the interpolation point coordinates.\n");

		// Zero out the interactionResult array.  Probably not necessary
		for (i = 0; i < numberOfTargets; i++) {
			interactionResult[i] = 0.0;
		    }


		// Fill the matrix of target - interpolation point kernel evaluations.  Note, this can/should be replaced with a threaded implementation on CPU or GPU.
        double dx, dy, dz, r;
//#pragma omp parallel for private(j,dx,dy,dz,r)
        for (i = 0; i < numberOfTargets; i++){

        	for (j = 0; j < numberOfInterpolationPoints; j++){

        		// Compute x, y, and z distances between target i and interpolation point j
        		dx = xT[ batch_ind[0] - 1 + i] - interpolationX[j];
        		dy = yT[ batch_ind[0] - 1 + i] - interpolationY[j];
        		dz = zT[ batch_ind[0] - 1 + i] - interpolationZ[j];
        		r = sqrt( dx*dx + dy*dy + dz*dz);
//        		if (r<1e-14){ printf("r is very small for a cluster-approx.  Why?\n");}
        		// Evaluate Kernel, store in kernelMatrix[i][j]
        		kernelMatrix[i*numberOfInterpolationPoints + j] = 1.0 / r;

        		// Perform the singularity subtraction piece here
//        		EnP[batch_ind[0] - 1 + i] -= qT[ batch_ind[0] - 1 + i]*Weights2[j]*exp(-r*r/kappaSq) / r;
        		EnP[batch_ind[0] - 1 + i] +=  ( Weights1[j] - Weights2[j]*qT[ batch_ind[0] - 1 + i]*exp(-r*r/kappaSq) ) / r;

//        		if (EnP[batch_ind[0] - 1 + i] > 1e6){printf("EnP blew up for target %d.  r = %12.5f, qT = %12.5f, Weight2 = %12.5f.\n",batch_ind[0] - 1 + i, r, qT[ batch_ind[0] - 1 + i],Weights2[j]);}

        	}

        }


//        // Multiply with CBLAS
//        double alpha=1;
//        double beta=0;
//
//        int incX = 1;
//        int incY = 1;
//
////        printf("Calling CBLAS_DGEMV.\n");
//        // Call CBAS for the f_j*w_j piece
//        cblas_dgemv(CblasRowMajor, CblasNoTrans, numberOfTargets, numberOfInterpolationPoints,
//        		alpha, kernelMatrix, numberOfInterpolationPoints, Weights1, incX, beta, interactionResult, incY );
//
//
//
//        // Add result to EnP, starting at index batch_ind[0] - 1
////		printf("Batch starting at: %d\n", batch_ind[0]-1);
//		for (i = 0; i < numberOfTargets; i++){
////			printf("Interation Result entry %d: %12.5e\n", batch_ind[0] - 1 + i, interactionResult[i]);
//			EnP[batch_ind[0] - 1 + i] += interactionResult[i];
//		}

		mkl_free(kernelMatrix);
		mkl_free(interactionResult);
		mkl_free(Weights1);
		mkl_free(Weights2);
		free(interpolationX);
		free(interpolationY);
		free(interpolationZ);


        
    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */


        if (p->num_children == 0) {
//        	printf("MAC rejected, and node has no children.  Calling pc_comp_dierct()...\n");
            pc_comp_direct_coulomb_SS(p->ibeg, p->iend, batch_ind[0], batch_ind[1],
                           xS, yS, zS, qS, wS, xT, yT, zT, qT, kappaSq, EnP);
        } else {
//        	printf("MAC rejected, recursing over children...\n");
            for (i = 0; i < p->num_children; i++) {
            	compute_pc_coulomb_SS(p->child[i], batch_ind, batch_mid, batch_rad,
                			xS, yS, zS, qS, wS, xT, yT, zT, qT, kappaSq, EnP);
            }
        }
    }

    return;

} /* END of function compute_pc */




/*
 * comp_direct directly computes the potential on the targets in the current
 * cluster due to the current source, determined by the global variable TARPOS
 */
void pc_comp_direct_coulomb_SS(int ibeg, int iend, int batch_ibeg, int batch_iend,
                    double *xS, double *yS, double *zS, double *qS, double *wS,
                    double *xT, double *yT, double *zT, double *qT, double kappaSq, double *EnP)
{
    /* local variables */
    int i, ii;
    double tx, ty, tz;

    double d_peng, r;

//#pragma omp parallel for private(i, d_peng,tx,ty,tz,r)
    for (ii = batch_ibeg - 1; ii < batch_iend; ii++) {
        d_peng = 0.0;
        for (i = ibeg - 1; i < iend; i++) {
            tx = xS[i] - xT[ii];
            ty = yS[i] - yT[ii];
            tz = zS[i] - zT[ii];
            r = sqrt(tx*tx + ty*ty + tz*tz);
//    		if (r<1e-10){ printf("r is very small for a direct sum piece, target %d and source %d\n",ii,i);}

            if (r > 1e-10){

            	d_peng += wS[i]* ( qS[i] - qT[ii]* exp(-r*r/kappaSq) )  / r;
//                V_Hartree_new[globalID] += weight_s * (rho_s -   rho_t * exp(- r*r / alphasq )   ) / r  # increment the new wavefunction value

//                if (d_peng > 1e10){printf("d_peng just blew up for target %d when interacting with source %d. r = %12.5f\n",ii,i,r);}

            }
        }
//        printf("d_peng from direct sum: %12.5e\n", d_peng);
        EnP[ii] += d_peng;
//        if (EnP[ii] > 1e5){printf("EnP just blew up for target %d when interacting with source %d\n",ii,i);}
    }

    return;

} /* END function pc_comp_direct_yuk */

