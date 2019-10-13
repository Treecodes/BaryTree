/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"



void pc_treecode_yuk_SS(struct tnode *p, struct batch *batches,
                     struct particles *sources, struct particles *targets, struct particles *clusters,
                     double kappa, double *tpeng, double *EnP)
{
    /* local variables */
    int i, j;
    
    for (i = 0; i < targets->num; i++){
        EnP[i] = 4.0*M_PI*targets->q[i]/kappa/kappa;  // 4*pi*f_t/k**2
    }




		double *EnP2;
		make_vector(EnP2,targets->num);
		for (i = 0; i < targets->num; i++)
			EnP2[i] = 0.0;

#pragma acc data copyin(sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], sources->w[0:sources->num], \
		targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num], targets->q[0:targets->num], \
		clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], clusters->q[0:clusters->num], clusters->w[0:clusters->num]) \
		copy(EnP2[0:targets->num])
    {
    for (i = 0; i < batches->num; i++) {
        for (j = 0; j < p->num_children; j++) {
            compute_pc_yuk_SS(p->child[j],
                batches->index[i], batches->center[i], batches->radius[i],
                sources->x, sources->y, sources->z, sources->q, sources->w,
                targets->x, targets->y, targets->z, targets->q, kappa, EnP2,
				clusters->x, clusters->y, clusters->z, clusters->q, clusters->w);
        }
    }
	#pragma acc wait
    } // end acc data region
    for (int k = 0; k < targets->num; k++){
    		if (EnP2[k] != 0.0){
    			EnP[k] += EnP2[k];
    		}
    		}
    free_vector(EnP2);
    

    *tpeng = sum(EnP, targets->num);
    
    return;
    
} /* END of function pc_treecode */




void compute_pc_yuk_SS(struct tnode *p,
                int *batch_ind, double *batch_mid, double batch_rad,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT, double kappa, double *EnP,
				double * clusterX, double * clusterY, double * clusterZ, double * clusterM, double * clusterM2 )
{
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

	int smallEnoughLeaf=0;
	if (p->numpar < torderlim*torderlim*torderlim){
		smallEnoughLeaf=1;
	}else{
		smallEnoughLeaf=0;
	}
	if (((p->radius + batch_rad) < dist * sqrt(thetasq)) && (p->sqradius != 0.00) && (smallEnoughLeaf==0)  ) {


		int numberOfTargets = batch_ind[1] - batch_ind[0] + 1;
		int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
		int clusterStart = numberOfInterpolationPoints*p->node_index;



		double xi,yi,zi,qi,r;
		int batchStart = batch_ind[0] - 1;

		int streamID = rand() % 2;
		# pragma acc kernels async(streamID) present(xT,yT,zT,qT,EnP, clusterX, clusterY, clusterZ, clusterM, clusterM2)
		{
		#pragma acc loop independent
		for (i = 0; i < numberOfTargets; i++){
			tempPotential = 0.0;
			xi = xT[ batchStart + i];
			yi = yT[ batchStart + i];
			zi = zT[ batchStart + i];
			qi = qT[ batchStart + i];
			#pragma acc loop independent
			for (j = 0; j < numberOfInterpolationPoints; j++){

				// Compute x, y, and z distances between target i and interpolation point j
				dxt = xi - clusterX[clusterStart + j];
				dyt = yi - clusterY[clusterStart + j];
				dzt = zi - clusterZ[clusterStart + j];

				r = sqrt(dxt*dxt + dyt*dyt + dzt*dzt);

//				tempPotential += (clusterM[clusterStart + j]-qi ) * exp(-kappa*r) / r;
				tempPotential += (clusterM[clusterStart + j]-qi*clusterM2[clusterStart + j] ) * exp(-kappa*r) / r;
	//			tempPotential += clusterM[clusterStart + j]* exp(-kappa*r) / r - qi*clusterM2[clusterStart + j] * exp(-kappa*r) / r;

							}
			#pragma acc atomic
			EnP[batchStart + i] += tempPotential;
		}
		}

        
    } else {
    /*
     * If MAC fails check to see if there are children. If not, perform direct
     * calculation. If there are children, call routine recursively for each.
     */


        if ( (p->num_children == 0)|(smallEnoughLeaf==1) ) {
//        	printf("MAC rejected, and node has no children.  Calling pc_comp_dierct()...\n");
            pc_comp_direct_yuk_SS(p->ibeg, p->iend, batch_ind[0], batch_ind[1],
                           xS, yS, zS, qS, wS, xT, yT, zT, qT, kappa, EnP);
        } else {
//        	printf("MAC rejected, recursing over children...\n");
            for (i = 0; i < p->num_children; i++) {
            	compute_pc_yuk_SS(p->child[i], batch_ind, batch_mid, batch_rad,
                			xS, yS, zS, qS, wS, xT, yT, zT, qT, kappa, EnP,
							clusterX, clusterY, clusterZ, clusterM, clusterM2);
            }
        }
    }

    return;

} /* END of function compute_pc */




/*
 * comp_direct directly computes the potential on the targets in the current
 * cluster due to the current source, determined by the global variable TARPOS
 */
void pc_comp_direct_yuk_SS(int ibeg, int iend, int batch_ibeg, int batch_iend,
                    double *xS, double *yS, double *zS, double *qS, double *wS,
                    double *xT, double *yT, double *zT, double *qT, double kappa, double *EnP)
{
    /* local variables */
    int i, ii;
    double tx, ty, tz;
    double d_peng, r;



    int streamID = rand() % 2;
	# pragma acc kernels async(streamID) present(xS,yS,zS,qS,wS,xT,yT,zT,qT,EnP)
    {
	#pragma acc loop independent
    for (ii = batch_ibeg - 1; ii < batch_iend; ii++) {
        d_peng = 0.0;
		#pragma acc loop independent
        for (i = ibeg - 1; i < iend; i++) {
            tx = xS[i] - xT[ii];
            ty = yS[i] - yT[ii];
            tz = zS[i] - zT[ii];
            r = sqrt(tx*tx + ty*ty + tz*tz);
            if (r > DBL_MIN) {
            	d_peng += (qS[i] - qT[ii]) * wS[i] * exp(-kappa*r) / r;
            }
        }
		#pragma acc atomic
        EnP[ii] += d_peng;
    }
    }

    return;

} /* END function pc_comp_direct_yuk */



void pc_interaction_list_treecode_yuk_SS(struct tnode_array *tree_array, struct particles *clusters, struct batch *batches,
                                  int *tree_inter_list, int *direct_inter_list,
                                  struct particles *sources, struct particles *targets,
                                  double *tpeng, double kappa, double *EnP)
{
        int i, j;
        int rank; int numProcs;	int ierr;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        int tree_numnodes = tree_array->numnodes;

//        for (i = 0; i < targets->num; i++)
//            EnP[i] = 0.0;

        for (i = 0; i < targets->num; i++){
			EnP[i] = 4.0*M_PI*targets->q[i]/kappa/kappa;  // 4*pi*f_t/k**2
		}



            double *EnP2, *EnP3;
            make_vector(EnP2,targets->num);
            make_vector(EnP3,targets->num);

            for (i = 0; i < targets->num; i++) {
                EnP3[i] = 0.0;
                EnP2[i] = 0.0;
            }

            double *xS = sources->x;
            double *yS = sources->y;
            double *zS = sources->z;
            double *qS = sources->q;
            double *wS = sources->w;

            double *xT = targets->x;
            double *yT = targets->y;
            double *zT = targets->z;
            double *qT = targets->q;

            double *xC = clusters->x;
            double *yC = clusters->y;
            double *zC = clusters->z;
            double *qC = clusters->q;
            double *wC = clusters->w;

//          printf("\n\nInside compute region, clusters->q[0] = %f\n\n",clusters->q[0]);
//          printf("\n\nInside compute region, clusters->q[213599] = %f\n\n",clusters->q[213599]);

            int * ibegs = tree_array->ibeg;
            int * iends = tree_array->iend;

            int * clusterInd = tree_array->cluster_ind;

#ifdef OPENACC_ENABLED
        #pragma acc data copyin(xS[0:sources->num], yS[0:sources->num], zS[0:sources->num], \
                            qS[0:sources->num], wS[0:sources->num], \
                            xT[0:targets->num], yT[0:targets->num], zT[0:targets->num], qT[0:targets->num], \
                            xC[0:clusters->num], yC[0:clusters->num], zC[0:clusters->num], qC[0:clusters->num],wC[0:clusters->num], \
                            tree_inter_list[0:tree_numnodes*batches->num], direct_inter_list[0:batches->num * numleaves], \
                            ibegs[0:tree_numnodes], iends[0:tree_numnodes]) copy(EnP3[0:targets->num], EnP2[0:targets->num])
#endif
        {

        int batch_ibeg, batch_iend, node_index;
        double dist;
        double tx, ty, tz;
        int i, j, k, ii, jj;
        double dxt,dyt,dzt,tempPotential;
        double temp_i[torderlim], temp_j[torderlim], temp_k[torderlim];

        int source_start;
        int source_end;

        double d_peng, r;
        double xi,yi,zi,qi;

        int numberOfTargets;
        int numberOfInterpolationPoints = torderlim*torderlim*torderlim;
        int clusterStart, batchStart;

        int numberOfClusterApproximations, numberOfDirectSums;
        int streamID;

        for (i = 0; i < batches->num; i++) {
            batch_ibeg = batches->index[i][0];
            batch_iend = batches->index[i][1];
            numberOfClusterApproximations = batches->index[i][2];
            numberOfDirectSums = batches->index[i][3];

//            printf("Rank %i, batch %i, number of cluster approximations: %i\n", rank, i, numberOfClusterApproximations);
//            printf("Rank %i, batch %i, number of direct interactions: %i\n", rank, i, numberOfDirectSums);

            numberOfTargets = batch_iend - batch_ibeg + 1;
            batchStart =  batch_ibeg - 1;

            for (j = 0; j < numberOfClusterApproximations; j++) {
                node_index = tree_inter_list[i * tree_numnodes + j];
//                clusterStart = numberOfInterpolationPoints*node_index;
                clusterStart = numberOfInterpolationPoints*clusterInd[node_index];

                streamID = j%3;
#ifdef OPENACC_ENABLED
                #pragma acc kernels async(streamID) //present(xT,yT,zT,qT,EnP, clusterX, clusterY, clusterZ, clusterM)
                {
                #pragma acc loop independent
#endif
                for (ii = 0; ii < numberOfTargets; ii++) {
                    tempPotential = 0.0;
                    xi = xT[ batchStart + ii];
                    yi = yT[ batchStart + ii];
                    zi = zT[ batchStart + ii];
                    qi = qT[ batchStart + ii];;

                    for (jj = 0; jj < numberOfInterpolationPoints; jj++) {
                        // Compute x, y, and z distances between target i and interpolation point j
                        dxt = xi - xC[clusterStart + jj];
                        dyt = yi - yC[clusterStart + jj];
                        dzt = zi - zC[clusterStart + jj];
                        r=sqrt(dxt*dxt + dyt*dyt + dzt*dzt);
//                        tempPotential += qC[clusterStart + jj]*exp(-kappa*r) / r;
                        tempPotential += ( qC[clusterStart + jj]-qi*wC[clusterStart + jj] )*exp(-kappa*r) / r;
//                        (qS[i] - qT[ii])

                    }
#ifdef OPENACC_ENABLED
                    #pragma acc atomic
#endif
                    EnP3[batchStart + ii] += tempPotential;
                }
#ifdef OPENACC_ENABLED
                } // end kernel
#endif
            } // end for loop over cluster approximations

            for (j = 0; j < numberOfDirectSums; j++) {

                node_index = direct_inter_list[i * tree_numnodes + j];

                source_start=ibegs[node_index]-1;
                source_end=iends[node_index];

                streamID = j%3;
#ifdef OPENACC_ENABLED
                # pragma acc kernels async(streamID)
                {
                #pragma acc loop independent
#endif
                for (ii = batchStart; ii < batchStart+numberOfTargets; ii++) {
                    d_peng = 0.0;

                    for (jj = source_start; jj < source_end; jj++) {
                        tx = xS[jj] - xT[ii];
                        ty = yS[jj] - yT[ii];
                        tz = zS[jj] - zT[ii];
                        r = sqrt(tx*tx + ty*ty + tz*tz);

                        if (r > DBL_MIN) {
                            d_peng += ( qS[jj] - qT[ii] ) * wS[jj]*exp(-kappa*r) / r;
                        }
                    }
#ifdef OPENACC_ENABLED
                    #pragma acc atomic
#endif
                    EnP2[ii] += d_peng;
                }
#ifdef OPENACC_ENABLED
                } // end kernel
#endif
            } // end loop over number of leaves
        } // end loop over target batches
#ifdef OPENACC_ENABLED
        #pragma acc wait
#endif
        } // end acc data region

        double totalDueToApprox = 0.0, totalDueToDirect = 0.0;
        totalDueToApprox = sum(EnP3, targets->num);
        totalDueToDirect = sum(EnP2, targets->num);
//        printf("Potential due to approximations: %f\n",totalDueToApprox);
//        printf("Potential due to direct: %f\n",totalDueToDirect);
        for (int k = 0; k < targets->num; k++) {
            if (EnP2[k] != 0.0)
                EnP[k] += EnP2[k];
                EnP[k] += EnP3[k];
            }

            free_vector(EnP2);
            free_vector(EnP3);

        *tpeng = sum(EnP, targets->num);

        return;

    } /* END of function pc_treecode */



