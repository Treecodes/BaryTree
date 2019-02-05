/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "batch.h"
#include "particles.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"

#include "mkl.h"
#include <omp.h>


void pc_create_tree_n0(struct tnode **p, struct particles *sources,
                       int ibeg, int iend, int maxparnode, double *xyzmm,
                       int level)
{

    /*local variables*/
    double x_mid, y_mid, z_mid, xl, yl, zl, lmax, t1, t2, t3;
    int i, j, loclev, numposchild;
    
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


    /* increment number of nodes */
    (*p)->node_index = numnodes;
    numnodes++;

    /* set node fields: number of particles, exist_ms, and xyz bounds */
    (*p)->numpar = iend - ibeg + 1;
    (*p)->exist_ms = 0;
    
    (*p)->x_min = minval(sources->x + ibeg - 1, (*p)->numpar);
    (*p)->x_max = maxval(sources->x + ibeg - 1, (*p)->numpar);
    (*p)->y_min = minval(sources->y + ibeg - 1, (*p)->numpar);
    (*p)->y_max = maxval(sources->y + ibeg - 1, (*p)->numpar);
    (*p)->z_min = minval(sources->z + ibeg - 1, (*p)->numpar);
    (*p)->z_max = maxval(sources->z + ibeg - 1, (*p)->numpar);
    
    //(*p)->x_min = xyzmm[0];
    //(*p)->x_max = xyzmm[1];
    //(*p)->y_min = xyzmm[2];
    //(*p)->y_max = xyzmm[3];
    //(*p)->z_min = xyzmm[4];
    //(*p)->z_max = xyzmm[5];
    

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


    if (maxlevel < level) maxlevel = level;

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

        pc_partition_8(sources->x, sources->y, sources->z, sources->q, sources->w,
                       xyzmms, xl, yl, zl, lmax, &numposchild,
                       x_mid, y_mid, z_mid, ind);

        loclev = level + 1;

        for (i = 0; i < numposchild; i++) {
            if (ind[i][0] <= ind[i][1]) {
                (*p)->num_children = (*p)->num_children + 1;

                for (j = 0; j < 6; j++)
                    lxyzmm[j] = xyzmms[j][i];
                
                pc_create_tree_n0(&((*p)->child[(*p)->num_children - 1]),
                                  sources, ind[i][0], ind[i][1],
                                  maxparnode, lxyzmm, loclev);
            }
        }
        
    } else {
        
        if (level < minlevel) minlevel = level;
        if (minpars > (*p)->numpar) minpars = (*p)->numpar;
        if (maxpars < (*p)->numpar) maxpars = (*p)->numpar;
        
        /* increment number of leaves */
        numleaves++;
    }

    return;

} /* END of function create_tree_n0 */



void pc_create_tree_array(struct tnode *p, struct tnode_array *tree_array)
{
    int i;

    /*midpoint coordinates, RADIUS and SQRADIUS*/
    tree_array->x_mid[p->node_index] = p->x_mid;
    tree_array->y_mid[p->node_index] = p->y_mid;
    tree_array->z_mid[p->node_index] = p->z_mid;

    tree_array->ibeg[p->node_index] = p->ibeg;
    tree_array->iend[p->node_index] = p->iend;

    for (i = 0; i < p->num_children; i++) {
        pc_create_tree_array(p->child[i], tree_array);
    }

    return;

} /* END of function create_tree_n0 */



void pc_partition_8(double *x, double *y, double *z, double *q, double *w, double xyzmms[6][8],
                    double xl, double yl, double zl, double lmax, int *numposchild,
                    double x_mid, double y_mid, double z_mid, int ind[8][2])
{
    /* local variables */
    int temp_ind, i, j;
    double critlen;

    *numposchild = 1;
    critlen = lmax / sqrt(2.0);

    if (xl >= critlen) {

        pc_partition(x, y, z, q, w, orderarr, ind[0][0], ind[0][1],
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
            pc_partition(y, x, z, q, w, orderarr, ind[i][0], ind[i][1],
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
            pc_partition(z, x, y, q, w, orderarr, ind[i][0], ind[i][1],
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

} /* END of function partition_8 */



void fill_in_cluster_data(struct particles *clusters, struct particles *sources, struct tnode *troot, int order){

	int pointsPerCluster = (order+1)*(order+1)*(order+1);
	int numInterpPoints = numnodes * pointsPerCluster;
	make_vector(clusters->x, numInterpPoints);
	make_vector(clusters->y, numInterpPoints);
	make_vector(clusters->z, numInterpPoints);
	make_vector(clusters->q, numInterpPoints);
	clusters->num=numInterpPoints;


	addNodeToArray( clusters->x, clusters->y, clusters->z, clusters->q, troot, sources, order, numInterpPoints, pointsPerCluster);
	return;
}

void addNodeToArray(double *x, double *y, double *z, double *q, struct tnode *p, struct particles *sources, int order, int numInterpPoints, int pointsPerCluster)
{
	int torderlim = order+1;
	int startingIndex = p->node_index * pointsPerCluster;
	int i;

	// compute moments
	make_vector(p->ms, (torderlim)*(torderlim)*(torderlim));
	make_vector(p->ms2, (torderlim)*(torderlim)*(torderlim));
	make_vector(p->tx, torderlim);
	make_vector(p->ty, torderlim);
	make_vector(p->tz, torderlim);

	for (i = 0; i < (torderlim)*(torderlim)*(torderlim); i++){
		p->ms[i] = 0.0;
		p->ms2[i] = 0.0;
	}

	if (p->node_index>8){ // experiment with not computing moments for top of tree, which are unlikely to be needed.
	pc_comp_ms(p, sources->x, sources->y, sources->z, sources->q, sources->w);
		p->exist_ms = 1;


		// fill in arrays, starting at startingIndex
		int k1,k2,k3;
		int kk = -1;
		for (k3 = 0; k3 < torderlim; k3++) {
			for (k2 = 0; k2 < torderlim; k2++) {
				for (k1 = 0; k1 < torderlim; k1++) {
					kk++;
					q[kk+startingIndex] = p->ms[kk];
					x[kk+startingIndex] = p->tx[k1];
					y[kk+startingIndex] = p->ty[k2];
					z[kk+startingIndex] = p->tz[k3];
				}
			}
		}
	}
	for (i = 0; i < p->num_children; i++) {
		addNodeToArray(x,y,z,q,p->child[i],sources,order,numInterpPoints,pointsPerCluster);
	}

	return;
}


void pc_treecode(struct tnode *p, struct batch *batches,
                 struct particles *sources, struct particles *targets, struct particles *clusters,
                 double *tpeng, double *EnP)
{
//	printf("Entered pc_treecoode.\n");
    /* local variables */
    int i, j;

    for (i = 0; i < targets->num; i++)
        EnP[i] = 0.0;
    
#pragma acc data region copyin(sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num], \
		targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num], targets->q[0:targets->num], EnP[0:targets->num], \
		clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num]) \
		copyout(EnP[0:targets->num])
    {

    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < p->num_children; j++) {
            compute_pc(p->child[j],
                batches->index[i], batches->center[i], batches->radius[i],
                sources->x, sources->y, sources->z, sources->q, sources->w,
                targets->x, targets->y, targets->z, targets->q, EnP,
				clusters->x, clusters->y, clusters->z, clusters->q);
        }
    }
}

    *tpeng = sum(EnP, targets->num);

    return;

} /* END of function pc_treecode */




void compute_pc(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double *EnP,
				double * clusterX, double * clusterY, double * clusterZ, double * clusterM )
{
//	printf("Entering compute_cp.  Batch start: %d.  Batch end: %d.\n", batch_ind[0]-1, batch_ind[1]);
    /* local variables */
    double dist;
    double tx, ty, tz;
    int i, j, k, ii, kk;
    double dxt,dyt,dzt,tempPotential;
    double temp_i[torderlim], temp_j[torderlim], temp_k[torderlim];

    /* determine DIST for MAC test */
    tx = batch_mid[0] - p->x_mid;
    ty = batch_mid[1] - p->y_mid;
    tz = batch_mid[2] - p->z_mid;
    dist = sqrt(tx*tx + ty*ty + tz*tz);


    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00) && (torderlim*torderlim*torderlim < p->numpar) ) {
//	if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00) ) {
    /*
     * If MAC is accepted and there is more than n0 particles
     * in the box, use the expansion for the approximation.
     */
//    	printf("MAC accepted, performing particle-cluster approximation...\n");

//        if (p->exist_ms == 0) {
//            make_vector(p->ms, (torderlim)*(torderlim)*(torderlim));
//            make_vector(p->ms2, (torderlim)*(torderlim)*(torderlim));
//            make_vector(p->tx, torderlim);
//            make_vector(p->ty, torderlim);
//            make_vector(p->tz, torderlim);
//
////            printf("Allocated vectors for ms, tx, ty, tz.\n");
//
//
//            for (i = 0; i < (torderlim)*(torderlim)*(torderlim); i++){
//                p->ms[i] = 0.0;
//            	p->ms2[i] = 0.0;
//            }
////            printf("Zeroed out p->ms \n");
//
//            pc_comp_ms(p, xS, yS, zS, qS, wS);
////            pc_comp_weights(p);
//            p->exist_ms = 1;
////            printf("Created 'moments' for node.\n");
//        }


	int numberOfTargets = batch_ind[1] - batch_ind[0] + 1;
	int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
	int clusterStart = numberOfInterpolationPoints*p->node_index;

//	double clusterX[numberOfInterpolationPoints], clusterY[numberOfInterpolationPoints], clusterZ[numberOfInterpolationPoints], localMoments[numberOfInterpolationPoints];
//	double clusterXYZM[4*numberOfInterpolationPoints];


//	// Fill some local cluster arrays from the cluster itself
//	int k1,k2,k3;
//	kk = -1;
//	for (k3 = 0; k3 < torderlim; k3++) {
//		for (k2 = 0; k2 < torderlim; k2++) {
//			for (k1 = 0; k1 < torderlim; k1++) {
//				kk++;
//				localMoments[kk] = p->ms[kk];
//				clusterX[kk] = p->tx[k1];
//				clusterY[kk] = p->ty[k2];
//				clusterZ[kk] = p->tz[k3];
//			}
//		}
//	}





	double xi,yi,zi;
	int batchStart = batch_ind[0] - 1;
//	double pointvals[4];
# pragma acc region present(xT,yT,zT,qT,EnP, clusterX, clusterY, clusterZ, clusterM)
    {
	#pragma acc loop independent
	for (i = 0; i < numberOfTargets; i++){
		tempPotential = 0.0;
		xi = xT[ batchStart + i];
		yi = yT[ batchStart + i];
		zi = zT[ batchStart + i];
		for (j = 0; j < numberOfInterpolationPoints; j++){

			// Compute x, y, and z distances between target i and interpolation point j
			dxt = xi - clusterX[clusterStart + j];
			dyt = yi - clusterY[clusterStart + j];
			dzt = zi - clusterZ[clusterStart + j];

//			pointvals = clusterXYZM[4*j];
//			dxt = xi - pointvals[0];
//			dyt = yi - pointvals[1];
//			dzt = zi - pointvals[2];
//			dxt = xi - clusterXYZM[4*j+0];
//			dyt = yi - clusterXYZM[4*j+1];
//			dzt = zi - clusterXYZM[4*j+2];


//			tempPotential += localMoments[j] / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);
			tempPotential += clusterM[clusterStart + j] / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);

						}

		EnP[batchStart + i] += tempPotential;
	}
    }
//        for (ii = batch_ind[0] - 1; ii < batch_ind[1]; ii++) {
//            tempPotential = 0.0;
//            kk = -1;
//            for (k = 0; k < torderlim; k++) {
//            	dzt = zT[ii] - clusterX[k];
//                for (j = 0; j < torderlim; j++) {
//                	dyt = yT[ii] - clusterY[j];
//                    for (i = 0; i < torderlim; i++) {
//                    	dxt = xT[ii] - clusterZ[i];
//                        kk++;
////                        tempPotential += localMoments[kk] / sqrt(temp_i[i] + temp_j[j] + temp_k[k]);
//                        tempPotential += localMoments[kk] / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);
////                        tempPotential += localMoments[k*torderlim*torderlim + j*torderlim + i] / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);
//
//                    }
//                }
//            }
//            EnP[ii] += tempPotential;
//        }
//    }


//        for (ii = batch_ind[0] - 1; ii < batch_ind[1]; ii++) {
//
//            for (i = 0; i < torderlim; i++) {
//                dxt = xT[ii] - p->tx[i];
//                dyt = yT[ii] - p->ty[i];
//                dzt = zT[ii] - p->tz[i];
//                temp_i[i] = dxt * dxt;
//                temp_j[i] = dyt * dyt;
//                temp_k[i] = dzt * dzt;
//            }
//
//            kk = -1;
//            for (k = 0; k < torderlim; k++) {
//                for (j = 0; j < torderlim; j++) {
//                    for (i = 0; i < torderlim; i++) {
//                        kk++;
//                        EnP[ii] += p->ms[kk] / sqrt(temp_i[i] + temp_j[j] + temp_k[k]);
//                    }
//                }
//            }
//        }



        // Allocate the matrices and vectors
//        printf("Working on batch from %d to %d\n\n", batch_ind[0] - 1,batch_ind[1]);
//        int numberOfTargets = batch_ind[1] - batch_ind[0] + 1;
//        int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
//
//        double *kernelMatrix 		= (double *)mkl_malloc(numberOfTargets * numberOfInterpolationPoints * sizeof(double),64);
////        double *kernelMatrix 		= (double *)malloc(numberOfTargets * numberOfInterpolationPoints * sizeof(double));
//        double *interactionResult 	= (double *)mkl_malloc(numberOfTargets * sizeof(double),64);
////        double *interactionResult 	= (double *)malloc(numberOfTargets * sizeof(double));
//
//
//        double *interpolationX = (double *)malloc(numberOfInterpolationPoints * sizeof(double));
//        double *interpolationY = (double *)malloc(numberOfInterpolationPoints * sizeof(double));
//        double *interpolationZ = (double *)malloc(numberOfInterpolationPoints * sizeof(double));
//        double *Weights 	   = (double *)mkl_malloc(numberOfInterpolationPoints * sizeof(double),64);
//        double *Weights 	   = (double *)malloc(numberOfInterpolationPoints * sizeof(double));

//        if (kernelMatrix == NULL || interactionResult == NULL || Weights == NULL) {
//			printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
//			mkl_free(kernelMatrix);
//			mkl_free(interactionResult);
//			mkl_free(Weights);
//			free(interpolationX);
//			free(interpolationY);
//			free(interpolationZ);
//        }
//
//
//        // Fill in the interpolation point coordinate vectors.  Not necessary, but helps modularize the next step, filling the kernel matrix.
//
//        int kk = -1;
//        int k1, k2, k3;
//		for (k3 = 0; k3 < torderlim; k3++) {
//			for (k2 = 0; k2 < torderlim; k2++) {
//				for (k1 = 0; k1 < torderlim; k1++) {
//					kk++;
//					interpolationX[kk] = p->tx[k1];
//					interpolationY[kk] = p->ty[k2];
//					interpolationZ[kk] = p->tz[k3];
//					Weights[kk] = p->ms[kk];
//
//				}
//			}
//		}
////		printf("Filled in the interpolation point coordinates.\n");
//
//		// Zero out the interactionResult array.  Probably not necessary
//		for (i = 0; i < numberOfTargets; i++) {
//			interactionResult[i] = 0.0;
//		    }
//
//
//		// Fill the matrix of target - interpolation point kernel evaluations.  Note, this can/should be replaced with a threaded implementation on CPU or GPU.
//        double dx, dy, dz, xi,yi,zi;
//#pragma omp parallel for private(j,dx,dy,dz,xi,yi,zi)
//        for (i = 0; i < numberOfTargets; i++){
//        	xi = xT[ batch_ind[0] - 1 + i];
//        	yi = yT[ batch_ind[0] - 1 + i];
//        	zi = zT[ batch_ind[0] - 1 + i];
//        	for (j = 0; j < numberOfInterpolationPoints; j++){
//
//        		// Compute x, y, and z distances between target i and interpolation point j
//        		dx = xi - interpolationX[j];
//        		dy = yi - interpolationY[j];
//        		dz = zi - interpolationZ[j];
//
//        		// Evaluate Kernel, store in kernelMatrix[i][j]
//        		kernelMatrix[i*numberOfInterpolationPoints + j] = 1 / sqrt( dx*dx + dy*dy + dz*dz);
//
//        	}
//
//        }
//		printf("Filled in the interaction matrix.\n");

        // Multiply kernel matrix with the vector of cluster weights.  Note, this can/should be replaced with a BLAS or cuBLAS call.
//        double tempSum;
//
//#pragma omp parallel for private(j,tempSum)
//        for (i = 0; i < numberOfTargets; i++){
//
//        	tempSum = 0.0;
//
//			for (j = 0; j < numberOfInterpolationPoints; j++){
//
//				tempSum += kernelMatrix[i*numberOfInterpolationPoints + j] * Weights[j];
//			}
//
//			interactionResult[i] = tempSum;
//
//        }
//
//
////        // Multiply with CBLAS
//        double alpha=1;
//        double beta=0;
//        int incX = 1;
//        int incY = 1;
//
//        cblas_dgemv(CblasRowMajor, CblasNoTrans, numberOfTargets, numberOfInterpolationPoints,
//        		alpha, kernelMatrix, numberOfInterpolationPoints, Weights, incX, beta, interactionResult, incY );
//
//
//
//        // Add result to EnP, starting at index batch_ind[0] - 1
////		printf("Batch starting at: %d\n", batch_ind[0]-1);
//		for (i = 0; i < numberOfTargets; i++){
////			printf("Interation Result entry %d: %12.5e\n", batch_ind[0] - 1 + i, interactionResult[i]);
//			EnP[batch_ind[0] - 1 + i] += interactionResult[i];
//		}
//
//		mkl_free(kernelMatrix);
//		mkl_free(interactionResult);
//		mkl_free(Weights);
//		free(interpolationX);
//		free(interpolationY);
//		free(interpolationZ);



    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */


        if (p->num_children == 0) {
//        	printf("MAC rejected, and node has no children.  Calling pc_comp_dierct()...\n");
            pc_comp_direct(p->ibeg, p->iend, batch_ind[0], batch_ind[1],
                           xS, yS, zS, qS, wS, xT, yT, zT, qT, EnP);
        } else {
//        	printf("MAC rejected, recursing over children...\n");
            for (i = 0; i < p->num_children; i++) {
                compute_pc(p->child[i], batch_ind, batch_mid, batch_rad,
                			xS, yS, zS, qS, wS, xT, yT, zT, qT, EnP,
							clusterX, clusterY, clusterZ, clusterM);
            }
        }
    }

    return;

} /* END of function compute_pc */

//void compute_pc(struct tnode *p,
//                int *batch_ind, double *batch_mid, double batch_rad,
//                double *xS, double *yS, double *zS, double *qS, double *wS,
//                double *xT, double *yT, double *zT, double *qT, double *EnP)
//{
//    /* local variables */
//    double dxt, dyt, dzt, dist;
//    double tx, ty, tz;
//    int i, j, k, kk, ii;
//    double temp_i[torderlim], temp_j[torderlim], temp_k[torderlim];
//
//    /* determine DIST for MAC test */
//    tx = batch_mid[0] - p->x_mid;
//    ty = batch_mid[1] - p->y_mid;
//    tz = batch_mid[2] - p->z_mid;
//    dist = sqrt(tx*tx + ty*ty + tz*tz);
//
//    int torder3 = torderlim*torderlim*torderlim;
//
//    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00)
//       && (torder3 < p->numpar)) {
//    /*
//     * If MAC is accepted and there is more than n0 particles
//     * in the box, use the expansion for the approximation.
//     */
//
//        if (p->exist_ms == 0) {
//            make_vector(p->ms, torder3);
//            make_vector(p->tx, torderlim);
//            make_vector(p->ty, torderlim);
//            make_vector(p->tz, torderlim);
//
//
//            for (i = 0; i < torder3; i++)
//                p->ms[i] = 0.0;
//
//            pc_comp_ms(p, xS, yS, zS, qS, wS);
//            p->exist_ms = 1;
//        }
//
//        for (ii = batch_ind[0] - 1; ii < batch_ind[1]; ii++) {
//
//            for (i = 0; i < torderlim; i++) {
//                dxt = xT[ii] - p->tx[i];
//                dyt = yT[ii] - p->ty[i];
//                dzt = zT[ii] - p->tz[i];
//                temp_i[i] = dxt * dxt;
//                temp_j[i] = dyt * dyt;
//                temp_k[i] = dzt * dzt;
//            }
//
//            kk = -1;
//            for (k = 0; k < torderlim; k++) {
//                for (j = 0; j < torderlim; j++) {
//                    for (i = 0; i < torderlim; i++) {
//                        kk++;
//                        EnP[ii] += p->ms[kk] / sqrt(temp_i[i] + temp_j[j] + temp_k[k]);
//                    }
//                }
//            }
//        }
//
//    } else {
//    /*
//     * If MAC fails check to see if there are children. If not, perform direct
//     * calculation. If there are children, call routine recursively for each.
//     */
//        if ((p->num_children == 0)) {
//            pc_comp_direct(p->ibeg, p->iend, batch_ind[0], batch_ind[1],
//                           xS, yS, zS, qS, wS, xT, yT, zT, qT, EnP);
//        } else {
//            for (i = 0; i < p->num_children; i++) {
//                compute_pc(p->child[i], batch_ind, batch_mid, batch_rad,
//                		xS, yS, zS, qS, wS, xT, yT, zT, qT, EnP);
//            }
//        }
//    }
//
//    return;
//
//} /* END of function compute_pc */




/*
 * comp_direct directly computes the potential on the targets in the current
 * cluster due to the current source, determined by the global variable TARPOS
 */
//void pc_comp_direct(int ibeg, int iend, int batch_ibeg, int batch_iend,
//                    __restrict__ double *xS,__restrict__ double *yS, __restrict__ double *zS, __restrict__ double *qS, __restrict__ double *wS,
//					__restrict__ double *xT, __restrict__ double *yT, __restrict__ double *zT, __restrict__ double *qT, __restrict__ double *EnP)
void pc_comp_direct(int ibeg, int iend, int batch_ibeg, int batch_iend,
                     double *xS, double *yS,  double *zS,  double *qS,  double *wS,
					 double *xT,  double *yT,  double *zT,  double *qT,  double *EnP)
{
    /* local variables */
    int i, ii;
    double tx, ty, tz;

    int batch_start=batch_ibeg - 1;
    int batch_end = batch_iend;

    int source_start=ibeg - 1;
    int source_end=iend;

    double d_peng, r;
# pragma acc region present(xS,yS,zS,qS,wS,xT,yT,zT,qT,EnP)
    {
	#pragma acc loop independent
    for (ii = batch_start; ii < batch_end; ii++) {
        d_peng = 0.0;
        for (i = source_start; i < source_end; i++) {
            tx = xS[i] - xT[ii];
            ty = yS[i] - yT[ii];
            tz = zS[i] - zT[ii];
            r = sqrt(tx*tx + ty*ty + tz*tz);
            if (r > DBL_MIN) {
            	d_peng += qS[i] * wS[i] / r;
            }
        }
        EnP[ii] += d_peng;
    }
    }
    return;

} /* END function pc_comp_direct */




//void pc_comp_direct(int ibeg, int iend, int batch_ibeg, int batch_iend,
//                        double *xS, double *yS, double *zS, double *qS, double *wS,
//                        double *xT, double *yT, double *zT, double *qT, double *EnP)
//    {
//
//
//    // Matrix-vector product for direct sum calculation, same as particle-cluster approx.
//	int i,j;
//    int numberOfTargets = batch_iend - batch_ibeg +1;  // need this +1?
//	int numberOfClusterPoints = iend - ibeg +1;
//
//	double *kernelMatrix 		= (double *)mkl_malloc(numberOfTargets * numberOfClusterPoints * sizeof(double),64);
//	double *interactionResult 	= (double *)mkl_malloc(numberOfTargets * sizeof(double),64);
//
//
//	double *clusterX = (double *)malloc(numberOfClusterPoints * sizeof(double));
//	double *clusterY = (double *)malloc(numberOfClusterPoints * sizeof(double));
//	double *clusterZ = (double *)malloc(numberOfClusterPoints * sizeof(double));
//	double *Weights	 = (double *)mkl_malloc(numberOfClusterPoints * sizeof(double),64);
//
//	// Zero out the interactionResult array.  Probably not necessary
//	for (i = 0; i < numberOfTargets; i++) {
//		interactionResult[i] = 0.0;
//		}
//
//	// Fill in cluster arrays so xS, yS, etc. only need to be accessed once per entry
//	for (j = 0; j < numberOfClusterPoints; j++) {
//		clusterX[j] = xS[ibeg-1+j];
//		clusterY[j] = yS[ibeg-1+j];
//		clusterZ[j] = zS[ibeg-1+j];
//		Weights[j] = qS[ibeg-1+j]*wS[ibeg-1+j];
//		}
//
//
//	// Fill the matrix of target - interpolation point kernel evaluations.  Note, this can/should be replaced with a threaded implementation on CPU or GPU.
//	double dx, dy, dz, xi, yi, zi, r;
//#pragma omp parallel for private(xi,yi,zi,j,dx,dy,dz,r)
//	for (i = 0; i < numberOfTargets; i++){
//		xi = xT[ batch_ibeg-1 + i]; // check these -1's
//		yi = yT[ batch_ibeg-1 + i];
//		zi = zT[ batch_ibeg-1 + i];
//		for (j = 0; j < numberOfClusterPoints; j++){
//
//			// Compute x, y, and z distances between target i and interpolation point j
//			dx = xi - clusterX[j];
//			dy = yi - clusterY[j];
//			dz = zi - clusterZ[j];
//
//			r = sqrt( dx*dx + dy*dy + dz*dz);
//			if (r > DBL_MIN) {
//			// Evaluate Kernel, store in kernelMatrix[i][j]
//				kernelMatrix[i*numberOfClusterPoints + j] = 1 / r;
//			}else{
//				kernelMatrix[i*numberOfClusterPoints + j] = 0.0;
//			}
//
//		}
//
//	}
//
//
//	// Multiply with CBLAS
//	double alpha=1;
//	double beta=0;
//	int incX = 1;
//	int incY = 1;
//
//	cblas_dgemv(CblasRowMajor, CblasNoTrans, numberOfTargets, numberOfClusterPoints,
//			alpha, kernelMatrix, numberOfClusterPoints, Weights, incX, beta, interactionResult, incY );
//
//	for (i = 0; i < numberOfTargets; i++){
//			EnP[batch_ibeg-1 + i] += interactionResult[i];
//	}
//
//	mkl_free(kernelMatrix);
//	mkl_free(interactionResult);
//	mkl_free(Weights);
//	free(clusterX);
//	free(clusterY);
//	free(clusterZ);
//
//
//    return;
//
//} /* END function pc_comp_direct */




/*
 * cp_comp_ms computes the moments for node p needed in the Taylor approximation
 */
void pc_comp_ms(struct tnode *p, double *x, double *y, double *z, double *q, double *w)
{

    int i, j, k1, k2, k3, kk;
    int a1exactind, a2exactind, a3exactind;
//    double dx, dy, dz, tx, ty, tz, qloc;
    double x0, x1, y0, y1, z0, z1;
    double sumA1, sumA2, sumA3;
    double xx, yy, zz;
    double *xibeg, *yibeg, *zibeg, *qibeg, *wibeg;
    
    double w1i[torderlim], w2j[torderlim], w3k[torderlim], dj[torderlim];
    double *Dd; //, **a1i, **a2j, **a3k;
    double *a1i, *a2j, *a3k;
    double *node_x, *node_y, *node_z;
    

    xibeg = &(x[p->ibeg-1]);
    yibeg = &(y[p->ibeg-1]);
    zibeg = &(z[p->ibeg-1]);
    qibeg = &(q[p->ibeg-1]);
    wibeg = &(w[p->ibeg-1]);
    
    x0 = p->x_min;
    x1 = p->x_max;
    y0 = p->y_min;
    y1 = p->y_max;
    z0 = p->z_min;
    z1 = p->z_max;
    
    for (i = 0; i < torderlim; i++) {
        p->tx[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);  // what is tt?
        p->ty[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        p->tz[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);

//        p->wx[i] = 0.0;
//        p->wy[i] = 0.0;
//        p->wz[i] = 0.0;  // the product wx[i]*wy[j]*wz[k] will give the quadrature weight at interpolation point (i,j,k)

    }
    
//    make_matrix(a1i, torderlim, p->numpar);
//    make_matrix(a2j, torderlim, p->numpar);
//    make_matrix(a3k, torderlim, p->numpar);

    make_vector(a1i, torderlim *p->numpar);
    make_vector(a2j, torderlim * p->numpar);
	make_vector(a3k, torderlim * p->numpar);

    make_vector(Dd, p->numpar);
    make_vector(node_x, torderlim);
    make_vector(node_y, torderlim);
    make_vector(node_z, torderlim);
    
    dj[0] = 0.5;
    dj[torder] = 0.5;
    for (j = 1; j < torder; j++)
        dj[j] = 1.0;
    
    for (j = 0; j < torderlim; j++) {
        w1i[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
        w2j[j] = w1i[j];
        w3k[j] = w1i[j];
    }
    
    sumA1 = 0.0;
    sumA2 = 0.0;
    sumA3 = 0.0;
    
    int pointsInNode = p->numpar;
    for (i=0;i<torderlim;i++){
    	node_x[i] = p->tx[i];
    	node_y[i] = p->ty[i];
    	node_z[i] = p->tz[i];
    }

//#pragma acc parallel loop independent
    for (i = 0; i < pointsInNode; i++) {
        xx = xibeg[i];
        yy = yibeg[i];
        zz = zibeg[i];
        
        a1exactind = -1;
        a2exactind = -1;
        a3exactind = -1;
//#pragma acc loop independent
        for (j = 0; j < torderlim; j++) {
//            a1i[j][i] = w1i[j] / (xx - node_x[j]);
            a1i[i*pointsInNode + j] = w1i[j] / (xx - node_x[j]);
//            a2j[j][i] = w2j[j] / (yy - node_y[j]);
            a2j[i*pointsInNode + j] = w2j[j] / (yy - node_y[j]);
//            a3k[j][i] = w3k[j] / (zz - node_z[j]);
            a3k[i*pointsInNode + j] = w3k[j] / (zz - node_z[j]);
            
            //if (isinf(a2j[j][i - p->ibeg + 1]))
            //    printf("a2j %d, %d: %f, %f, %f\n", i, j, a2j[j][i - p->ibeg + 1], yy, p->ty[j]);
            
            //if (isinf(a1i[j][i - p->ibeg + 1]))
            //    printf("a1i %d, %d: %f, %f, %f\n", i, j, a1i[j][i - p->ibeg + 1], xx, p->tx[j]);

            //if (isinf(a3k[j][i - p->ibeg + 1]))
            //    printf("a3k %d, %d: %f, %f, %f\n", i, j, a3k[j][i - p->ibeg + 1], zz, p->tz[j]);

//            sumA1 += a1i[j][i];
//            sumA2 += a2j[j][i];
//            sumA3 += a3k[j][i];

            sumA1 += a1i[i*pointsInNode + j];
			sumA2 += a2j[i*pointsInNode + j];
			sumA3 += a3k[i*pointsInNode + j];
            
            if (fabs(xx - node_x[j]) < DBL_MIN) a1exactind = j;
            if (fabs(yy - node_y[j]) < DBL_MIN) a2exactind = j;
            if (fabs(zz - node_z[j]) < DBL_MIN) a3exactind = j;
        }
        
        if (a1exactind > -1) {
            sumA1 = 1.0;
            for (j = 0; j < torderlim; j++) a1i[i*pointsInNode + j] = 0.0;
            a1i[i*pointsInNode + a1exactind] = 1.0;
        }
        
        if (a2exactind > -1) {
            sumA2 = 1.0;
            for (j = 0; j < torderlim; j++) a2j[i*pointsInNode + j] = 0.0;
            a2j[i*pointsInNode + a2exactind] = 1.0;
        }
        
        if (a3exactind > -1) {
            sumA3 = 1.0;
            for (j = 0; j < torderlim; j++) a3k[i*pointsInNode + j] = 0.0;
            a3k[i*pointsInNode + a3exactind] = 1.0;
        }
        
        Dd[i] = 1.0 / (sumA1 * sumA2 * sumA3);
        
        sumA1 = 0.0;
        sumA2 = 0.0;
        sumA3 = 0.0;
    }
    
//    kk = -1;

//#pragma acc parallel loop
    for (k3 = 0; k3 < torderlim; k3++) {
        for (k2 = 0; k2 < torderlim; k2++) {
            for (k1 = 0; k1 < torderlim; k1++) {
//                kk++;
                for (i = 0; i < p->numpar; i++) {
//                	d_ms[k3*torderlim*torderlim + k2*torderlim + k1] += a1i[k1][i] * a2j[k2][i] * a3k[k3][i] * Dd[i] * qibeg[i] * wibeg[i] ;  // in this case, multiple the function value by the quadrature weight.  Revisit for sing. subt.
                	p->ms[k3*torderlim*torderlim + k2*torderlim + k1] += a1i[i*pointsInNode + k1] * a2j[i*pointsInNode + k2] * a3k[i*pointsInNode + k3] * Dd[i] * qibeg[i] * wibeg[i] ;  // in this case, multiple the function value by the quadrature weight.  Revisit for sing. subt.
//                	d_ms2[k3*torderlim*torderlim + k2*torderlim + k1] += a1i[k1][i] * a2j[k2][i] * a3k[k3][i] * Dd[i] *  wibeg[i] ;  // Anterpolate only the quadrature weights w_i for ms2, instead of w_i*f_i
                    p->ms2[k3*torderlim*torderlim + k2*torderlim + k1] += a1i[i*pointsInNode + k1] * a2j[i*pointsInNode + k2] * a3k[i*pointsInNode + k3] * Dd[i] *  wibeg[i] ;  // Anterpolate only the quadrature weights w_i for ms2, instead of w_i*f_i
                    
                    //if (p->ms[kk] != p->ms[kk])
                    //    printf("%d, %d, %d, %d: %f, %f, %f, %f, %f, %f\n", k1, k2, k3, i,
                    //    p->ms[kk], Dd[i], q[p->ibeg-1+i], a1i[k1][i], a2j[k2][i], a3k[k3][i]);
                }
            }
        }
    }
    

    free_vector(a1i);
    free_vector(a2j);
    free_vector(a3k);
    free_vector(Dd);
    
    return;
    
} /* END function cp_comp_ms */






void pc_make_interaction_list(struct tnode *p, struct batch *batches,
                              int **tree_inter_list, int **direct_inter_list)
{
    /* local variables */
    int i, j;
    int tree_index_counter;
    int direct_index_counter;

    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < numnodes; j++) {
            tree_inter_list[i][j] = -1;
        }
        for (j = 0; j < numleaves; j++) {
            direct_inter_list[i][j] = -1;
        }
    }
    
    for (i = 0; i < batches->num; i++) {
        tree_index_counter = 0;
        direct_index_counter = 0;
        
        pc_compute_interaction_list(p,
                batches->index[i], batches->center[i], batches->radius[i],
                tree_inter_list[i], direct_inter_list[i],
                &tree_index_counter, &direct_index_counter);
    }

    return;

} /* END of function pc_treecode */




void pc_compute_interaction_list(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                int *batch_tree_list, int *batch_direct_list,
                int *tree_index_counter, int *direct_index_counter)
{
    /* local variables */
    double tx, ty, tz, dist;
    int i;

    /* determine DIST for MAC test */
    tx = batch_mid[0] - p->x_mid;
    ty = batch_mid[1] - p->y_mid;
    tz = batch_mid[2] - p->z_mid;
    dist = sqrt(tx*tx + ty*ty + tz*tz);

    if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00)) {
    /*
     * If MAC is accepted and there is more than 1 particle
     * in the box, use the expansion for the approximation.
     */

        batch_tree_list[*tree_index_counter] = p->node_index;
        (*tree_index_counter)++;

    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */
        if (p->num_children == 0) {
            batch_direct_list[*direct_index_counter] = p->node_index;
            (*direct_index_counter)++;

        } else {
            for (i = 0; i < p->num_children; i++) {
                pc_compute_interaction_list(p->child[i], batch_ind, batch_mid, batch_rad,
                           batch_tree_list, batch_direct_list,
                           tree_index_counter, direct_index_counter);
            }
        }
    }

    return;

} /* END of function compute_pc */
