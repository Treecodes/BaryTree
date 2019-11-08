/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "array.h"
#include "globvars.h"
#include "tnode.h"
#include "batch.h"
#include "particles.h"
#include "tools.h"

#include "kernels/kernels.h"
#include "kernels/coulomb.h"
#include "kernels/yukawa.h"
#include "kernels/coulombSingularitySubtraction.h"
#include "kernels/yukawaSingularitySubtraction.h"
#include "tree.h"



void fill_in_cluster_data(struct particles *clusters, struct particles *sources, struct tnode *troot,
                          int interpolationOrder, struct tnode_array * tree_array)
{
    int tree_numnodes = tree_array->numnodes;
    int interpolationPointsPerCluster = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);
    int totalNumberInterpolationPoints = tree_numnodes * interpolationPointsPerCluster;

    make_vector(clusters->x, totalNumberInterpolationPoints);
    make_vector(clusters->y, totalNumberInterpolationPoints);
    make_vector(clusters->z, totalNumberInterpolationPoints);
    make_vector(clusters->q, totalNumberInterpolationPoints);
    make_vector(clusters->w, totalNumberInterpolationPoints);  // will be used in singularity subtraction
    clusters->num = totalNumberInterpolationPoints;

    for (int i = 0; i < totalNumberInterpolationPoints; i++) {
        clusters->x[i] = 0.0;
        clusters->y[i] = 0.0;
        clusters->z[i] = 0.0;
        clusters->q[i] = 0.0;
        clusters->w[i] = 1.0;
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
        double *wC = clusters->w;

        int totalNumberSourcePoints = sources->num;
        int interpolationPointsPerDimension = (interpolationOrder+1);


#ifdef OPENACC_ENABLED
        #pragma acc data copyin(tt[0:interpolationPointsPerDimension], \
        xS[0:totalNumberSourcePoints], yS[0:totalNumberSourcePoints], zS[0:totalNumberSourcePoints], qS[0:totalNumberSourcePoints], wS[0:totalNumberSourcePoints]) \
        copy(xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], zC[0:totalNumberInterpolationPoints], qC[0:totalNumberInterpolationPoints], wC[0:totalNumberInterpolationPoints] )
        {
#endif
            for (int i = 0; i < tree_numnodes; i++) {
                pc_comp_ms_modifiedF(tree_array, i, interpolationOrder, xS, yS, zS, qS, wS, xC, yC, zC, qC, wC);
            }
#ifdef OPENACC_ENABLED
            #pragma acc wait
        } // end ACC DATA REGION
#endif

    return;
}


void pc_comp_ms_modifiedF(struct tnode_array *tree_array, int idx, int interpolationOrder,
        double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW)
{
    int interpolationPointsPerCluster = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);
    int sourcePointsInCluster = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInSourcesArray = tree_array->ibeg[idx]-1;

    double x0, x1, y0, y1, z0, z1;  // bounding box

    double weights[(interpolationOrder+1)];
    double dj[(interpolationOrder+1)];
    double *modifiedF, *modifiedF2;
    make_vector(modifiedF,sourcePointsInCluster);
    make_vector(modifiedF2,sourcePointsInCluster);

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
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < sourcePointsInCluster; j++) {
        modifiedF[j] = qS[startingIndexInSourcesArray+j] * wS[startingIndexInSourcesArray+j];
//        modifiedF[j] = qS[startingIndexInSourcesArray+j];
//        modifiedF2[j] = wS[startingIndexInSourcesArray+j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
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
    for (int j = 0; j < interpolationOrder+1; j++){
        dj[j] = 1.0;
        if (j==0) dj[j] = 0.5;
        if (j==interpolationOrder) dj[j]=0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < (interpolationOrder+1); j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

    // Compute modified f values
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
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
//        modifiedF2[i] /= denominator;

    }

    // Compute moments for each interpolation point
    double numerator, xn, yn, zn, temp, temp2;
    int k1, k2, k3, kk;
    double w1,w2,w3;

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
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
        temp2 = 0.0;
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
//            temp2 += numerator*modifiedF2[i];
//            temp2 += modifiedF2[i];
        }
        clusterQ[startingIndexInClustersArray + j] += temp;
//        clusterW[startingIndexInClustersArray + j] += temp2;
//        clusterW[startingIndexInClustersArray + j] = 1.0;
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



void pc_interaction_list_treecode(struct tnode_array *tree_array, struct batch *batches,
                                  int *tree_inter_list, int *direct_inter_list,
                                  double *xS, double *yS, double *zS, double *qS, double *wS,
                                  double *xT, double *yT, double *zT, double *qT,
                                  double *xC, double *yC, double *zC, double *qC, double *wC,
                                  double *totalPotential, double *pointwisePotential, int interpolationOrder,
                                  int numSources, int numTargets, int numClusters,
                                  int batch_approx_offset, int batch_direct_offset,
                                  char *kernelName, double kernel_parameter,
                                  char *approximationName)
{
        int rank, numProcs, ierr;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

        int tree_numnodes = tree_array->numnodes;

        double *potentialDueToDirect, *potentialDueToApprox;
        make_vector(potentialDueToDirect, numTargets);
        make_vector(potentialDueToApprox, numTargets);

        for (int i = 0; i < numTargets; i++) {
            potentialDueToApprox[i] = 0.0;
            potentialDueToDirect[i] = 0.0;
        }

        int *ibegs = tree_array->ibeg;
        int *iends = tree_array->iend;
        int *clusterInd = tree_array->cluster_ind;

#ifdef OPENACC_ENABLED
        #pragma acc data copyin(xS[0:numSources], yS[0:numSources], zS[0:numSources], \
                            qS[0:numSources], wS[0:numSources], \
                            xT[0:numTargets], yT[0:numTargets], zT[0:numTargets], qT[0:numTargets], \
                            xC[0:numClusters], yC[0:numClusters], zC[0:numClusters], qC[0:numClusters], \
                            tree_inter_list[0:batch_approx_offset*batches->num], \
                            direct_inter_list[0:batch_direct_offset*batches->num], \
                            ibegs[0:tree_numnodes], iends[0:tree_numnodes]) \
                            copy(potentialDueToApprox[0:numTargets], potentialDueToDirect[0:numTargets])
#endif
        {

        int numberOfInterpolationPoints = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);

        for (int i = 0; i < batches->num; i++) {
            int batch_ibeg = batches->index[i][0];
            int batch_iend = batches->index[i][1];
            int numberOfClusterApproximations = batches->index[i][2];
            int numberOfDirectSums = batches->index[i][3];

            int numberOfTargets = batch_iend - batch_ibeg + 1;
            int batchStart =  batch_ibeg - 1;


/**********************************************************/
/************** POTENTIAL FROM APPROX *********************/
/**********************************************************/

            for (int j = 0; j < numberOfClusterApproximations; j++) {
                int node_index = tree_inter_list[i * batch_approx_offset + j];
                int clusterStart = numberOfInterpolationPoints*clusterInd[node_index];

                int streamID = j%3;


        /***********************************************/
        /***************** Coulomb *********************/
        /***********************************************/

                if (strcmp(kernelName, "coulomb") == 0) {

                    if (strcmp(approximationName, "Lagrange") == 0) {
                
                        coulombApproximationLagrange(   numberOfTargets, numberOfInterpolationPoints, batchStart, clusterStart,
                                                        xT, yT, zT,
                                                        xC, yC, zC, qC,
                                                        potentialDueToApprox, streamID);
                    }else{
                        printf("Invalid approximationName.\n");
                        return;
                    }

        /***********************************************/
        /***************** Yukawa **********************/
        /***********************************************/

                } else if (strcmp(kernelName, "yukawa") == 0) {

                    if (strcmp(approximationName, "Lagrange") == 0) {

                        yukawaApproximationLagrange(   numberOfTargets, numberOfInterpolationPoints, batchStart, clusterStart,
                                                    xT, yT, zT,
                                                    xC, yC, zC, qC,
                                                    kernel_parameter, potentialDueToApprox, streamID);
                    }else{
                        printf("Invalid approximationName.\n");
                        return;
                    }

        /***********************************************/
        /***** Coulomb with Singularity Subtraction ****/
        /***********************************************/

                } else if (strcmp(kernelName, "coulombSingularitySubtraction") == 0) {

                    if (strcmp(approximationName, "Lagrange") == 0) {
                    coulombSingularitySubtractionApproximationLagrange(  numberOfTargets, numberOfInterpolationPoints, batchStart, clusterStart,
                                                                        xT, yT, zT, qT,
                                                                        xC, yC, zC, qC, wC,
                                                                        kernel_parameter, potentialDueToApprox, streamID);
                    }else{
                        printf("Invalid approximationName.\n");
                        return;
                    }

        /***********************************************/
        /***** Yukawa with Singularity Subtraction *****/
        /***********************************************/

                } else if (strcmp(kernelName, "yukawaSingularitySubtraction") == 0) {

                    if (strcmp(approximationName, "Lagrange") == 0) {

                        yukawaSingularitySubtractionApproximationLagrange(  numberOfTargets, numberOfInterpolationPoints, batchStart, clusterStart,
                                                                            xT, yT, zT, qT,
                                                                            xC, yC, zC, qC, wC,
                                                                            kernel_parameter, potentialDueToApprox, streamID);
                    }else{
                        printf("Invalid approximationName.\n");
                        return;
                    }

                } // end kernel selection

            } // end loop over cluster approximations



/**********************************************************/
/************** POTENTIAL FROM DIRECT *********************/
/**********************************************************/

            for (int j = 0; j < numberOfDirectSums; j++) {

                int node_index = direct_inter_list[i * batch_direct_offset + j];
                int source_start = ibegs[node_index]-1;
                int source_end = iends[node_index];
                int number_of_sources_in_cluster = source_end-source_start;
                int streamID = j%3;

        /***********************************************/
        /***************** Coulomb *********************/
        /***********************************************/

                if (strcmp(kernelName, "coulomb") == 0) {

                    coulombDirect(  numberOfTargets, number_of_sources_in_cluster, batchStart, source_start,
                                    xT, yT, zT,
                                    xS, yS, zS, qS, wS,
                                    potentialDueToDirect, streamID);

        /***********************************************/
        /***************** Yukawa **********************/
        /***********************************************/

                } else if (strcmp(kernelName, "yukawa") == 0) {

                    yukawaDirect(  numberOfTargets, number_of_sources_in_cluster, batchStart, source_start,
                                    xT, yT, zT,
                                    xS, yS, zS, qS, wS,
                                    kernel_parameter, potentialDueToDirect, streamID);

        /***********************************************/
        /***** Coulomb with Singularity Subtraction ****/
        /***********************************************/

                } else if (strcmp(kernelName, "coulombSingularitySubtraction") == 0) {

                    coulombSingularitySubtractionDirect( numberOfTargets, number_of_sources_in_cluster, batchStart, source_start,
                                                        xT, yT, zT, qT,
                                                        xS, yS, zS, qS, wS,
                                                        kernel_parameter, potentialDueToDirect, streamID);

        /***********************************************/
        /***** Yukawa with Singularity Subtraction *****/
        /***********************************************/

                } else if (strcmp(kernelName, "yukawaSingularitySubtraction") == 0) {

                    yukawaSingularitySubtractionDirect( numberOfTargets, number_of_sources_in_cluster, batchStart, source_start,
                                                        xT, yT, zT, qT,
                                                        xS, yS, zS, qS, wS,
                                                        kernel_parameter, potentialDueToDirect, streamID);



                } // end kernel selection

            } // end loop over number of leaves

        } // end loop over target batches
#ifdef OPENACC_ENABLED
        #pragma acc wait
#endif
        } // end acc data region


        double totalDueToApprox = sum(potentialDueToApprox, numTargets);
        double totalDueToDirect = sum(potentialDueToDirect, numTargets);

//        printf("Total due to direct = %f\n", totalDueToDirect);
//        printf("Total due to approx = %f\n", totalDueToApprox);

        for (int k = 0; k < numTargets; k++) {
            pointwisePotential[k] += potentialDueToDirect[k];
            pointwisePotential[k] += potentialDueToApprox[k];
         }

        free_vector(potentialDueToDirect);
        free_vector(potentialDueToApprox);

        *totalPotential = sum(pointwisePotential, numTargets);

        return;

    } /* END of function pc_treecode */
