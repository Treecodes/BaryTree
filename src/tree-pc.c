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
#include "batch.h"
#include "particles.h"
#include "tools.h"

#include "partition.h"
#include "tree.h"



void fill_in_cluster_data(struct particles *clusters, struct particles *sources, struct tnode *troot, int interpolationOrder, struct tnode_array * tree_array){

	int tree_numnodes = tree_array->numnodes;
    int interpolationPointsPerCluster = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);
    int totalNumberInterpolationPoints = tree_numnodes * interpolationPointsPerCluster;
    make_vector(clusters->x, totalNumberInterpolationPoints);
    make_vector(clusters->y, totalNumberInterpolationPoints);
    make_vector(clusters->z, totalNumberInterpolationPoints);
    make_vector(clusters->q, totalNumberInterpolationPoints);
    make_vector(clusters->w, totalNumberInterpolationPoints);  // will be used in singularity subtraction
    clusters->num=totalNumberInterpolationPoints;

    for (int i = 0; i < totalNumberInterpolationPoints; i++) {
        clusters->x[i]=0.0;
        clusters->y[i]=0.0;
        clusters->z[i]=0.0;
        clusters->q[i]=0.0;
        clusters->w[i]=0.0;
    }


        double *xS = sources->x;
        double *yS = sources->y;
        double *zS = sources->z;
        double *qS = sources->q;
        double *wS = sources->w;

        double *xC = clusters->x;
        double *yC = clusters->y;
        double *zC = clusters->z;
        double *qC = clusters->q;

        int totalNumberSourcePoints = sources->num;
        int interpolationPointsPerDimension = (interpolationOrder+1);

	printf("xs, ys, zs, qs: %f,%f,%f,%f\n",xS[0],yS[0],zS[0],qS[0]);
#ifdef OPENACC_ENABLED
        #pragma acc data copyin(tt[0:interpolationPointsPerDimension], \
        xS[0:totalNumberSourcePoints], yS[0:totalNumberSourcePoints], zS[0:totalNumberSourcePoints], qS[0:totalNumberSourcePoints], wS[0:totalNumberSourcePoints]) \
        copy(xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], zC[0:totalNumberInterpolationPoints], qC[0:totalNumberInterpolationPoints] )
        {
#endif
#ifdef OPENMP_ENABLED
        #pragma omp target enter data map(to: tt[0:interpolationPointsPerDimension], xS[0:totalNumberSourcePoints], yS[0:totalNumberSourcePoints], zS[0:totalNumberSourcePoints], qS[0:totalNumberSourcePoints], wS[0:totalNumberSourcePoints], xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], zC[0:totalNumberInterpolationPoints], qC[0:totalNumberInterpolationPoints] )
        {
#endif
            for (int i = 0; i < tree_numnodes; i++) {
            	pc_comp_ms_modifiedF(tree_array, i, interpolationOrder, xS, yS, zS, qS, wS, xC, yC, zC, qC);
            printf("i = %i\n", i);   
		 }
#ifdef OPENACC_ENABLED
            #pragma acc wait
        } // end ACC DATA REGION
#endif
#ifdef OPENMP_ENABLED
	//#pragma omp taskwait
	#pragma omp target exit data map(from: xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], zC[0:totalNumberInterpolationPoints], qC[0:totalNumberInterpolationPoints] ) \
					map(delete: xS[0:totalNumberSourcePoints], yS[0:totalNumberSourcePoints], zS[0:totalNumberSourcePoints], qS[0:totalNumberSourcePoints], wS[0:totalNumberSourcePoints])
	}
#endif
	printf("xc, yc, zc, qc: %f,%f,%f,%f\n",xC[0],yC[0],zC[0],qC[0]);
	printf("xc, yc, zc, qc: %f,%f,%f,%f\n",xC[200],yC[200],zC[200],qC[200]);
	printf("xs, ys, zs, qs: %f,%f,%f,%f\n",xS[0],yS[0],zS[0],qS[0]);
    return;
}


void pc_comp_ms_modifiedF(struct tnode_array * tree_array, int idx, int interpolationOrder,
        double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ)
{
    int interpolationPointsPerCluster = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);
    int sourcePointsInCluster = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInSourcesArray = tree_array->ibeg[idx]-1;

    double x0, x1, y0, y1, z0, z1;  // bounding box

    double weights[(interpolationOrder+1)];
    double dj[(interpolationOrder+1)];
    double *modifiedF;
    make_vector(modifiedF,sourcePointsInCluster);

    double nodeX[(interpolationOrder+1)], nodeY[(interpolationOrder+1)], nodeZ[(interpolationOrder+1)];

    int *exactIndX, *exactIndY, *exactIndZ;
    make_vector(exactIndX, sourcePointsInCluster);
    make_vector(exactIndY, sourcePointsInCluster);
    make_vector(exactIndZ, sourcePointsInCluster);

    x0 = tree_array->x_min[idx];
    x1 = tree_array->x_max[idx];
    y0 = tree_array->y_min[idx];
    y1 = tree_array->y_max[idx];
    z0 = tree_array->z_min[idx];
    z1 = tree_array->z_max[idx];

    // Make and zero-out arrays to store denominator sums
    double sumX, sumY, sumZ;
    double sx,sy,sz,cx,cy,cz,denominator,w;

    int streamID = rand() % 4;
#ifdef OPENACC_ENABLED
		#pragma acc kernels async(streamID) present(xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ, clusterQ,tt) \
		create(modifiedF[0:sourcePointsInCluster],exactIndX[0:sourcePointsInCluster],exactIndY[0:sourcePointsInCluster],exactIndZ[0:sourcePointsInCluster], \
				nodeX[0:(interpolationOrder+1)],nodeY[0:(interpolationOrder+1)],nodeZ[0:(interpolationOrder+1)],weights[0:(interpolationOrder+1)],dj[0:(interpolationOrder+1)])
		{
		#pragma acc loop independent

#endif
#ifdef OPENMP_ENABLED
    	#pragma omp target enter data map( alloc: xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ, clusterQ,tt, modifiedF[0:sourcePointsInCluster],exactIndX[0:sourcePointsInCluster],exactIndY[0:sourcePointsInCluster],exactIndZ[0:sourcePointsInCluster], \
            nodeX[0:(interpolationOrder+1)],nodeY[0:(interpolationOrder+1)],nodeZ[0:(interpolationOrder+1)],weights[0:(interpolationOrder+1)],dj[0:(interpolationOrder+1)])
		#pragma omp target teams //nowait
    	{
#endif
#ifdef OPENMP_ENABLED
	#pragma omp distribute parallel for
#endif
    for (int j = 0; j < sourcePointsInCluster; j++) {
        modifiedF[j] = qS[startingIndexInSourcesArray+j] * wS[startingIndexInSourcesArray+j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
#ifdef OPENMP_ENABLED
        #pragma omp distribute parallel for
#endif
    for (int i = 0; i < (interpolationOrder+1); i++) {
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);

    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
#ifdef OPENMP_ENABLED
        #pragma omp distribute parallel for
#endif
    for (int j = 0; j < interpolationOrder+1; j++){
        dj[j] = 1.0;
        if (j==0) dj[j] = 0.5;
        if (j==interpolationOrder) dj[j]=0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
#ifdef OPENMP_ENABLED
        #pragma omp distribute parallel for
#endif
    for (int j = 0; j < (interpolationOrder+1); j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

    // Compute modified f values
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
#ifdef OPENMP_ENABLED
        #pragma omp distribute parallel for
#endif
    for (int i = 0; i < sourcePointsInCluster; i++) { // loop through the source points

        sumX=0.0;
        sumY=0.0;
        sumZ=0.0;

        sx = xS[startingIndexInSourcesArray+i];
        sy = yS[startingIndexInSourcesArray+i];
        sz = zS[startingIndexInSourcesArray+i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (int j = 0; j < (interpolationOrder+1); j++) {  // loop through the degree

            cx = sx-nodeX[j];
            cy = sy-nodeY[j];
            cz = sz-nodeZ[j];

            if (fabs(cx)<DBL_MIN) exactIndX[i]=j;
            if (fabs(cy)<DBL_MIN) exactIndY[i]=j;
            if (fabs(cz)<DBL_MIN) exactIndZ[i]=j;

            // Increment the sums
            w = weights[j];
            sumX += w / (cx);
            sumY += w / (cy);
            sumZ += w / (cz);

        }

        denominator = 1.0;
        if (exactIndX[i]==-1) denominator *= sumX;
        if (exactIndY[i]==-1) denominator *= sumY;
        if (exactIndZ[i]==-1) denominator *= sumZ;

        modifiedF[i] /= denominator;

    }

    // Compute moments for each interpolation point
    double numerator, xn, yn, zn, temp;
    int k1, k2, k3, kk;
    double w1,w2,w3;

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
#ifdef OPENMP_ENABLED
        #pragma omp distribute parallel for
#endif
    for (int j = 0; j < interpolationPointsPerCluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point
        // compute k1, k2, k3 from j
        k1 = j%(interpolationOrder+1);
        kk = (j-k1)/(interpolationOrder+1);
        k2 = kk%(interpolationOrder+1);
        kk = kk - k2;
        k3 = kk / (interpolationOrder+1);

        cz = nodeZ[k3];
        w3 = weights[k3];

        cy = nodeY[k2];
        w2 = weights[k2];

        cx = nodeX[k1];
        w1 = weights[k1];

        // Fill cluster X, Y, and Z arrays
        clusterX[startingIndexInClustersArray + j] = cx;
        clusterY[startingIndexInClustersArray + j] = cy;
        clusterZ[startingIndexInClustersArray + j] = cz;

        // Increment cluster Q array
        temp = 0.0;
#ifdef OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (int i = 0; i < sourcePointsInCluster; i++) {  // loop over source points
            sx = xS[startingIndexInSourcesArray+i];
            sy = yS[startingIndexInSourcesArray+i];
            sz = zS[startingIndexInSourcesArray+i];

            numerator=1.0;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
            if (exactIndX[i]==-1) {
                numerator *=  w1 / (sx - cx);
            } else {
                if (exactIndX[i]!=k1) numerator *= 0;
            }

            if (exactIndY[i]==-1) {
                numerator *=  w2 / (sy - cy);
            } else {
                if (exactIndY[i]!=k2) numerator *= 0;
            }

            if (exactIndZ[i]==-1) {
                numerator *=  w3 / (sz - cz);
            } else {
                if (exactIndZ[i]!=k3) numerator *= 0;
            }

            temp += numerator * modifiedF[i];
        }
        clusterQ[startingIndexInClustersArray + j] += temp;
    }
#ifdef OPENACC_ENABLED
    }
#endif
#ifdef OPENMP_ENABLED
    }
#pragma omp target exit data map( delete: modifiedF[0:sourcePointsInCluster],exactIndX[0:sourcePointsInCluster],exactIndY[0:sourcePointsInCluster],exactIndZ[0:sourcePointsInCluster], \
            nodeX[0:(interpolationOrder+1)],nodeY[0:(interpolationOrder+1)],nodeZ[0:(interpolationOrder+1)],weights[0:(interpolationOrder+1)],dj[0:(interpolationOrder+1)])
#endif

    free_vector(modifiedF);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}


//void pc_interaction_list_treecode(struct tnode_array *tree_array, struct particles *clusters, struct batch *batches,
//                                  int *tree_inter_list, int *direct_inter_list,
//                                  struct particles *sources, struct particles *targets,
//                                  double *totalPotential, double *pointwisePotential, int interpolationOrder)

void pc_interaction_list_treecode(struct tnode_array *tree_array, struct batch *batches,
                                  int *tree_inter_list, int *direct_inter_list,
								  double *xS, double *yS, double *zS, double *qS, double *wS,
								  double *xT, double *yT, double *zT, double *qT,
								  double *xC, double *yC, double *zC, double *qC,
                                  double *totalPotential, double *pointwisePotential, int interpolationOrder,
								  int numSources, int numTargets, int numClusters)
{
        int i, j;
        int rank; int numProcs;	int ierr;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        int tree_numnodes = tree_array->numnodes;

//        for (i = 0; i < targets->num; i++){
//            pointwisePotential[i] = 0.0;
//        }



            double *potentialDueToDirect, *potentialDueToApprox;
            make_vector(potentialDueToDirect,numTargets);
            make_vector(potentialDueToApprox,numTargets);

            for (i = 0; i < numTargets; i++) {
                potentialDueToApprox[i] = 0.0;
                potentialDueToDirect[i] = 0.0;
            }

//            double *xS = sources->x;
//            double *yS = sources->y;
//            double *zS = sources->z;
//            double *qS = sources->q;
//            double *wS = sources->w;
//
//            double *xT = targets->x;
//            double *yT = targets->y;
//            double *zT = targets->z;
//            double *qT = targets->q;
//
//            double *xC = clusters->x;
//            double *yC = clusters->y;
//            double *zC = clusters->z;
//            double *qC = clusters->q;

            int * ibegs = tree_array->ibeg;
            int * iends = tree_array->iend;

            int * clusterInd = tree_array->cluster_ind;

#ifdef OPENACC_ENABLED
        #pragma acc data copyin(xS[0:numSources], yS[0:numSources], zS[0:numSources], \
                            qS[0:numSources], wS[0:numSources], \
                            xT[0:numTargets], yT[0:numTargets], zT[0:numTargets], qT[0:numTargets], \
                            xC[0:numClusters], yC[0:numClusters], zC[0:numClusters], qC[0:numClusters], \
                            tree_inter_list[0:tree_numnodes*batches->num], direct_inter_list[0:batches->num * numleaves], \
                            ibegs[0:tree_numnodes], iends[0:tree_numnodes]) copy(potentialDueToApprox[0:numTargets], potentialDueToDirect[0:numTargets])
#endif
#ifdef OPENMP_ENABLED
        #pragma omp target enter data map(to: 	xS[0:numSources], yS[0:numSources], zS[0:numSources], qS[0:numSources], wS[0:numSources], \
												xT[0:numTargets], yT[0:numTargets], zT[0:numTargets], qT[0:numTargets], \
												xC[0:numClusters], yC[0:numClusters], zC[0:numClusters], qC[0:numClusters], \
												tree_inter_list[0:tree_numnodes*batches->num], direct_inter_list[0:batches->num * numleaves], \
												ibegs[0:tree_numnodes], iends[0:tree_numnodes], potentialDueToApprox[0:numTargets], potentialDueToDirect[0:numTargets])
#endif
        {

        int batch_ibeg, batch_iend, node_index;
        double dist;
        double tx, ty, tz;
        int i, j, k, ii, jj;
        double dxt,dyt,dzt,tempPotential;
        double temp_i[(interpolationOrder+1)], temp_j[(interpolationOrder+1)], temp_k[(interpolationOrder+1)];

        int source_start;
        int source_end;

        double d_peng, r;
        double xi,yi,zi;

        int numberOfTargets;
        int numberOfInterpolationPoints = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);
        int clusterStart, batchStart;

        int numberOfClusterApproximations, numberOfDirectSums;
        int streamID;

        for (i = 0; i < batches->num; i++) {
            batch_ibeg = batches->index[i][0];
            batch_iend = batches->index[i][1];
            numberOfClusterApproximations = batches->index[i][2];
            numberOfDirectSums = batches->index[i][3];

            numberOfTargets = batch_iend - batch_ibeg + 1;
            batchStart =  batch_ibeg - 1;

            for (j = 0; j < numberOfClusterApproximations; j++) {
                node_index = tree_inter_list[i * tree_numnodes + j];
                clusterStart = numberOfInterpolationPoints*clusterInd[node_index];

                streamID = j%3;
#ifdef OPENACC_ENABLED
                #pragma acc kernels async(streamID)
                {
                #pragma acc loop independent
#endif
#ifdef OPENMP_ENABLED
				#pragma omp target teams nowait
                #pragma omp  distribute parallel for
#endif
                for (ii = 0; ii < numberOfTargets; ii++) {
                    tempPotential = 0.0;
                    xi = xT[ batchStart + ii];
                    yi = yT[ batchStart + ii];
                    zi = zT[ batchStart + ii];

                    for (jj = 0; jj < numberOfInterpolationPoints; jj++) {
                        // Compute x, y, and z distances between target i and interpolation point j
                        dxt = xi - xC[clusterStart + jj];
                        dyt = yi - yC[clusterStart + jj];
                        dzt = zi - zC[clusterStart + jj];
                        tempPotential += qC[clusterStart + jj] / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);

                    }
#ifdef OPENACC_ENABLED
                    #pragma acc atomic // is this still needed now that we don't have openMP?  Or was this necessary due to streams?
#endif
#ifdef OPENMP_ENABLED
                    #pragma omp atomic // is this still needed now that we don't have openMP?  Or was this necessary due to streams?
#endif
                    potentialDueToApprox[batchStart + ii] += tempPotential;
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
#ifdef OPENMP_ENABLED
				#pragma omp target teams nowait
                #pragma omp  distribute parallel for
#endif
                for (ii = batchStart; ii < batchStart+numberOfTargets; ii++) {
                    d_peng = 0.0;

                    for (jj = source_start; jj < source_end; jj++) {
                        tx = xS[jj] - xT[ii];
                        ty = yS[jj] - yT[ii];
                        tz = zS[jj] - zT[ii];
                        r = sqrt(tx*tx + ty*ty + tz*tz);

                        if (r > DBL_MIN) {
                            d_peng += qS[jj] * wS[jj] / r;
                        }
                    }
#ifdef OPENACC_ENABLED
                    #pragma acc atomic
#endif
#ifdef OPENMP_ENABLED
                    #pragma omp atomic
#endif
                    potentialDueToDirect[ii] += d_peng;
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

#ifdef OPENMP_ENABLED
		#pragma omp target exit data map(from: potentialDueToDirect[0:numTargets], potentialDueToApprox[0:numTargets]) device(0)
#endif
        double totalDueToApprox = 0.0, totalDueToDirect = 0.0;
        totalDueToApprox = sum(potentialDueToApprox, numTargets);
        totalDueToDirect = sum(potentialDueToDirect, numTargets);
        //printf("Total due to direct = %f\n", totalDueToDirect);
        //printf("Total due to approx = %f\n", totalDueToApprox);
        for (int k = 0; k < numTargets; k++) {
//            if (potentialDueToDirect[k] != 0.0){
                pointwisePotential[k] += potentialDueToDirect[k];
                pointwisePotential[k] += potentialDueToApprox[k];
//        	}
         }

            free_vector(potentialDueToDirect);
            free_vector(potentialDueToApprox);

        *totalPotential = sum(pointwisePotential, numTargets);

        return;

    } /* END of function pc_treecode */
