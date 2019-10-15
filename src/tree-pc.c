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



void fill_in_cluster_data(struct particles *clusters, struct particles *sources, struct tnode *troot, int order, struct tnode_array * tree_array){

	int tree_numnodes = tree_array->numnodes;
    int interpolationPointsPerCluster = (order+1)*(order+1)*(order+1);
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


#ifdef OPENACC_ENABLED
        #pragma acc data copyin(tt[0:torderlim], \
        xS[0:totalNumberSourcePoints], yS[0:totalNumberSourcePoints], zS[0:totalNumberSourcePoints], qS[0:totalNumberSourcePoints], wS[0:totalNumberSourcePoints], \
        xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], zC[0:totalNumberInterpolationPoints], qZ[0:totalNumberInterpolationPoints])
        {
#endif
            for (int i = 0; i < tree_numnodes; i++) {

            	pc_comp_ms_modifiedF(tree_array, i, xS, yS, zS, qS, wS,
									 xC, yC, zC, qC);
            }
#ifdef OPENACC_ENABLED
            #pragma acc wait
        } // end ACC DATA REGION
#endif

    return;
}


void pc_comp_ms_modifiedF(struct tnode_array * tree_array, int idx,
        double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ)
{
    int i,j,k;
    int interpolationPointsPerCluster = torderlim*torderlim*torderlim;
    int pointsInNode = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClusters = idx * interpolationPointsPerCluster;
    int startingIndexInSources = tree_array->ibeg[idx]-1;

    double x0, x1, y0, y1, z0, z1;  // bounding box

    double weights[torderlim];
    double dj[torderlim];
    double *modifiedF;
    make_vector(modifiedF,pointsInNode);

    double nodeX[torderlim], nodeY[torderlim], nodeZ[torderlim];

    int *exactIndX, *exactIndY, *exactIndZ;
    make_vector(exactIndX, pointsInNode);
    make_vector(exactIndY, pointsInNode);
    make_vector(exactIndZ, pointsInNode);

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
    create(modifiedF[0:pointsInNode],exactIndX[0:pointsInNode],exactIndY[0:pointsInNode],exactIndZ[0:pointsInNode], \
            nodeX[0:torderlim],nodeY[0:torderlim],nodeZ[0:torderlim],weights[0:torderlim],dj[0:torderlim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (j = 0; j < pointsInNode; j++) {
        modifiedF[j] = qS[startingIndexInSources+j] * wS[startingIndexInSources+j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (i = 0; i < torderlim; i++) {
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);

    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (j = 0; j < torder+1; j++){
        dj[j] = 1.0;
        if (j==0) dj[j] = 0.5;
        if (j==torder) dj[j]=0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (j = 0; j < torderlim; j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

    // Compute modified f values
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (i = 0; i < pointsInNode; i++) { // loop through the source points

        sumX=0.0;
        sumY=0.0;
        sumZ=0.0;

        sx = xS[startingIndexInSources+i];
        sy = yS[startingIndexInSources+i];
        sz = zS[startingIndexInSources+i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (j = 0; j < torderlim; j++) {  // loop through the degree

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
    for (j = 0; j < interpolationPointsPerCluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point
        // compute k1, k2, k3 from j
        k1 = j%torderlim;
        kk = (j-k1)/torderlim;
        k2 = kk%torderlim;
        kk = kk - k2;
        k3 = kk / torderlim;

        cz = nodeZ[k3];
        w3 = weights[k3];

        cy = nodeY[k2];
        w2 = weights[k2];

        cx = nodeX[k1];
        w1 = weights[k1];

        // Fill cluster X, Y, and Z arrays
        clusterX[startingIndexInClusters + j] = cx;
        clusterY[startingIndexInClusters + j] = cy;
        clusterZ[startingIndexInClusters + j] = cz;

        // Increment cluster Q array
        temp = 0.0;
#ifdef OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (i = 0; i < pointsInNode; i++) {  // loop over source points
            sx = xS[startingIndexInSources+i];
            sy = yS[startingIndexInSources+i];
            sz = zS[startingIndexInSources+i];

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
        clusterQ[startingIndexInClusters + j] += temp;
    }
#ifdef OPENACC_ENABLED
    }
#endif

    free_vector(modifiedF);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}


void pc_interaction_list_treecode(struct tnode_array *tree_array, struct particles *clusters, struct batch *batches,
                                  int *tree_inter_list, int *direct_inter_list,
                                  struct particles *sources, struct particles *targets,
                                  double *tpeng, double *EnP)
{
        int i, j;
        int rank; int numProcs;	int ierr;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        int tree_numnodes = tree_array->numnodes;

        for (i = 0; i < targets->num; i++)
            EnP[i] = 0.0;



            double *potentialDueToDirect, *potentialDueToApprox;
            make_vector(potentialDueToDirect,targets->num);
            make_vector(potentialDueToApprox,targets->num);

            for (i = 0; i < targets->num; i++) {
                potentialDueToApprox[i] = 0.0;
                potentialDueToDirect[i] = 0.0;
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

//          printf("\n\nInside compute region, clusters->q[0] = %f\n\n",clusters->q[0]);
//          printf("\n\nInside compute region, clusters->q[213599] = %f\n\n",clusters->q[213599]);

            int * ibegs = tree_array->ibeg;
            int * iends = tree_array->iend;

            int * clusterInd = tree_array->cluster_ind;

#ifdef OPENACC_ENABLED
        #pragma acc data copyin(xS[0:sources->num], yS[0:sources->num], zS[0:sources->num], \
                            qS[0:sources->num], wS[0:sources->num], \
                            xT[0:targets->num], yT[0:targets->num], zT[0:targets->num], qT[0:targets->num], \
                            xC[0:clusters->num], yC[0:clusters->num], zC[0:clusters->num], qC[0:clusters->num], \
                            tree_inter_list[0:tree_numnodes*batches->num], direct_inter_list[0:batches->num * numleaves], \
                            ibegs[0:tree_numnodes], iends[0:tree_numnodes]) copy(potentialDueToApprox[0:targets->num], potentialDueToDirect[0:targets->num])
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
        double xi,yi,zi;

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

                    for (jj = 0; jj < numberOfInterpolationPoints; jj++) {
                        // Compute x, y, and z distances between target i and interpolation point j
                        dxt = xi - xC[clusterStart + jj];
                        dyt = yi - yC[clusterStart + jj];
                        dzt = zi - zC[clusterStart + jj];
                        tempPotential += qC[clusterStart + jj] / sqrt(dxt*dxt + dyt*dyt + dzt*dzt);

                    }
#ifdef OPENACC_ENABLED
                    #pragma acc atomic
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

        double totalDueToApprox = 0.0, totalDueToDirect = 0.0;
        totalDueToApprox = sum(potentialDueToApprox, targets->num);
        totalDueToDirect = sum(potentialDueToDirect, targets->num);
//        printf("Potential due to approximations: %f\n",totalDueToApprox);
//        printf("Potential due to direct: %f\n",totalDueToDirect);
        for (int k = 0; k < targets->num; k++) {
            if (potentialDueToDirect[k] != 0.0)
                EnP[k] += potentialDueToDirect[k];
                EnP[k] += potentialDueToApprox[k];
            }

            free_vector(potentialDueToDirect);
            free_vector(potentialDueToApprox);

        *tpeng = sum(EnP, targets->num);

        return;

    } /* END of function pc_treecode */
