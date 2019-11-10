
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

#include "clusters.h"

void pc_comp_ms_modifiedF(struct tnode_array *tree_array, int idx, int interpolationOrder,
                          double *xS, double *yS, double *zS, double *qS, double *wS,
                          double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW);

void pc_comp_ms_modifiedF_SS(struct tnode_array *tree_array, int idx, int interpolationOrder,
                          double *xS, double *yS, double *zS, double *qS, double *wS,
                          double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW);

void pc_comp_ms_modifiedF_hermite(struct tnode_array *tree_array, int idx, int interpolationOrder,
                          int totalNumberInterpolationPoints,
                          double *xS, double *yS, double *zS, double *qS, double *wS,
                          double *clusterX, double *clusterY, double *clusterZ,
                          double *clusterQ, double *clusterW);

void pc_comp_ms_modifiedF_hermite_SS(struct tnode_array *tree_array, int idx, int interpolationOrder,
                          int totalNumberInterpolationPoints,
                          double *xS, double *yS, double *zS, double *qS, double *wS,
                          double *clusterX, double *clusterY, double *clusterZ,
                          double *clusterQ, double *clusterW);



void Clusters_PC_Setup(struct particles *clusters, struct particles *sources,
                       int interpolationOrder, struct tnode_array *tree_array,
                       char *approxName, char *singularityHandling)
{
    int tree_numnodes = tree_array->numnodes;
    int totalNumberSourcePoints = sources->num;

    int interpOrderLim = interpolationOrder + 1;
    int interpolationPointsPerCluster = interpOrderLim * interpOrderLim * interpOrderLim;

    int totalNumberInterpolationPoints = tree_numnodes * interpolationPointsPerCluster;
    int totalNumberInterpolationCharges = totalNumberInterpolationPoints;
    int totalNumberInterpolationWeights = totalNumberInterpolationPoints; 

    clusters->num = totalNumberInterpolationPoints;
    make_vector(clusters->x, totalNumberInterpolationPoints);
    make_vector(clusters->y, totalNumberInterpolationPoints);
    make_vector(clusters->z, totalNumberInterpolationPoints);
    

    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->x[i] = 0.0;
    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->y[i] = 0.0;
    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->z[i] = 0.0;
    
    
    if (strcmp(approxName, "lagrange") == 0) {
        make_vector(clusters->q, totalNumberInterpolationCharges);
        make_vector(clusters->w, totalNumberInterpolationWeights);
        for (int i = 0; i < totalNumberInterpolationCharges; i++) clusters->q[i] = 0.0;
        
        if (strcmp(singularityHandling, "skipping") == 0) {
            for (int i = 0; i < totalNumberInterpolationWeights; i++) clusters->w[i] = 1.0;

        } else if (strcmp(singularityHandling, "subtraction") == 0) {
            for (int i = 0; i < totalNumberInterpolationWeights; i++) clusters->w[i] = 0.0;

        } else {
            exit(1);
        }
        
    } else if (strcmp(approxName, "hermite") == 0) {
        totalNumberInterpolationCharges *= 8;
        make_vector(clusters->q, totalNumberInterpolationCharges);
        for (int i = 0; i < totalNumberInterpolationCharges; i++) clusters->q[i] = 0.0;
        
        if (strcmp(singularityHandling, "skipping") == 0) {
            make_vector(clusters->w, totalNumberInterpolationWeights);
            for (int i = 0; i < totalNumberInterpolationWeights; i++) clusters->w[i] = 1.0;

        } else if (strcmp(singularityHandling, "subtraction") == 0) {
            totalNumberInterpolationWeights *= 8;
            make_vector(clusters->w, totalNumberInterpolationWeights);
            for (int i = 0; i < totalNumberInterpolationWeights; i++) clusters->w[i] = 0.0;

        } else {
            exit(1);
        }
        
    } else {
        exit(1);
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


#ifdef OPENACC_ENABLED
    #pragma acc data copyin(tt[0:interpOrderLim], ww[0:interpOrderLim], \
                            xS[0:totalNumberSourcePoints], yS[0:totalNumberSourcePoints], \
                            zS[0:totalNumberSourcePoints], qS[0:totalNumberSourcePoints], \
                            wS[0:totalNumberSourcePoints]) \
                       copy(xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], \
                            zC[0:totalNumberInterpolationPoints], qC[0:totalNumberInterpolationCharges], \
                            wC[0:totalNumberInterpolationWeights])
    {
#endif

    if ((strcmp(approxName, "lagrange") == 0) && (strcmp(singularityHandling, "skipping") == 0)) {
        for (int i = 0; i < tree_numnodes; i++)
            pc_comp_ms_modifiedF(tree_array, i, interpolationOrder, xS, yS, zS, qS, wS, xC, yC, zC, qC, wC);

    } else if ((strcmp(approxName, "lagrange") == 0) && (strcmp(singularityHandling, "subtraction") == 0)) {
        for (int i = 0; i < tree_numnodes; i++)
            pc_comp_ms_modifiedF_SS(tree_array, i, interpolationOrder, xS, yS, zS, qS, wS, xC, yC, zC, qC, wC);

    } else if ((strcmp(approxName, "hermite") == 0) && (strcmp(singularityHandling, "skipping") == 0)) {
        for (int i = 0; i < tree_numnodes; i++)
            pc_comp_ms_modifiedF_hermite(tree_array, i, interpolationOrder, totalNumberInterpolationPoints,
                                         xS, yS, zS, qS, wS, xC, yC, zC, qC, wC);

    } else if ((strcmp(approxName, "hermite") == 0) && (strcmp(singularityHandling, "subtraction") == 0)) {
        for (int i = 0; i < tree_numnodes; i++)
            pc_comp_ms_modifiedF_hermite_SS(tree_array, i, interpolationOrder, totalNumberInterpolationPoints,
                                            xS, yS, zS, qS, wS, xC, yC, zC, qC, wC);

    } else {
        exit(1);
    }
    
#ifdef OPENACC_ENABLED
    #pragma acc wait
    } // end ACC DATA REGION
#endif

    return;
}





/************************************/
/***** LOCAL FUNCTIONS **************/
/************************************/



void pc_comp_ms_modifiedF(struct tnode_array *tree_array, int idx, int interpolationOrder,
        double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW)
{

    int interpOrderLim = interpolationOrder + 1;
    int interpolationPointsPerCluster = interpOrderLim * interpOrderLim * interpOrderLim;
    int sourcePointsInCluster = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInSourcesArray = tree_array->ibeg[idx]-1;

    double weights[interpOrderLim], dj[interpOrderLim];
    double nodeX[interpOrderLim], nodeY[interpOrderLim], nodeZ[interpOrderLim];
    
    double *modifiedF, *modifiedF2;
    make_vector(modifiedF, sourcePointsInCluster);
    make_vector(modifiedF2, sourcePointsInCluster);

    int *exactIndX, *exactIndY, *exactIndZ;
    make_vector(exactIndX, sourcePointsInCluster);
    make_vector(exactIndY, sourcePointsInCluster);
    make_vector(exactIndZ, sourcePointsInCluster);


    double x0 = tree_array->x_min[idx];
    double x1 = tree_array->x_max[idx];
    double y0 = tree_array->y_min[idx];
    double y1 = tree_array->y_max[idx];
    double z0 = tree_array->z_min[idx];
    double z1 = tree_array->z_max[idx];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ, clusterQ, tt) \
                       create(modifiedF[0:sourcePointsInCluster], exactIndX[0:sourcePointsInCluster], \
                              exactIndY[0:sourcePointsInCluster], exactIndZ[0:sourcePointsInCluster], \
                              nodeX[0:(interpolationOrder+1)], nodeY[0:(interpolationOrder+1)], \
                              nodeZ[0:(interpolationOrder+1)], weights[0:(interpolationOrder+1)], \
                              dj[0:(interpolationOrder+1)])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < sourcePointsInCluster; j++) {
        modifiedF[j] = qS[startingIndexInSourcesArray + j] * wS[startingIndexInSourcesArray + j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interpOrderLim; i++) {
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interpOrderLim; j++) {
        dj[j] = 1.0;
        if (j == 0) dj[j] = 0.5;
        if (j == interpolationOrder) dj[j] = 0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interpOrderLim; j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

    // Compute modified f values
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < sourcePointsInCluster; i++) { // loop through the source points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double sx = xS[startingIndexInSourcesArray+i];
        double sy = yS[startingIndexInSourcesArray+i];
        double sz = zS[startingIndexInSourcesArray+i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:sumX) reduction(+:sumY) reduction(+:sumZ)
#endif
        for (int j = 0; j < (interpolationOrder+1); j++) {  // loop through the degree

            double cx = sx - nodeX[j];
            double cy = sy - nodeY[j];
            double cz = sz - nodeZ[j];

            if (fabs(cx) < DBL_MIN) exactIndX[i] = j;
            if (fabs(cy) < DBL_MIN) exactIndY[i] = j;
            if (fabs(cz) < DBL_MIN) exactIndZ[i] = j;

            // Increment the sums
            double w = weights[j];
            sumX += w / cx;
            sumY += w / cy;
            sumZ += w / cz;

        }

        double denominator = 1.0;
        if (exactIndX[i] == -1) denominator *= sumX;
        if (exactIndY[i] == -1) denominator *= sumY;
        if (exactIndZ[i] == -1) denominator *= sumZ;

        modifiedF[i] /= denominator;
    }


#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interpolationPointsPerCluster; j++) {
        int k1 = j%(interpolationOrder+1);
        int kk = (j-k1)/(interpolationOrder+1);
        int k2 = kk%(interpolationOrder+1);
        kk = kk - k2;
        int k3 = kk / (interpolationOrder+1);

        double cz = nodeZ[k3];
        double w3 = weights[k3];

        double cy = nodeY[k2];
        double w2 = weights[k2];

        double cx = nodeX[k1];
        double w1 = weights[k1];

        // Fill cluster X, Y, and Z arrays
        clusterX[startingIndexInClustersArray + j] = cx;
        clusterY[startingIndexInClustersArray + j] = cy;
        clusterZ[startingIndexInClustersArray + j] = cz;

        // Increment cluster Q array
        double temp = 0.0;
#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp)
#endif
        for (int i = 0; i < sourcePointsInCluster; i++) {  // loop over source points
            double sx = xS[startingIndexInSourcesArray + i];
            double sy = yS[startingIndexInSourcesArray + i];
            double sz = zS[startingIndexInSourcesArray + i];

            double numerator = 1.0;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
            if (exactIndX[i] == -1) {
                numerator *= w1 / (sx - cx);
            } else {
                if (exactIndX[i] != k1) numerator *= 0;
            }

            if (exactIndY[i] == -1) {
                numerator *= w2 / (sy - cy);
            } else {
                if (exactIndY[i] != k2) numerator *= 0;
            }

            if (exactIndZ[i] == -1) {
                numerator *= w3 / (sz - cz);
            } else {
                if (exactIndZ[i] != k3) numerator *= 0;
            }

            temp += numerator * modifiedF[i];

        }
        
        clusterQ[startingIndexInClustersArray + j] += temp;

    }
#ifdef OPENACC_ENABLED
    } //end acc kernels region
#endif

    free_vector(modifiedF);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}


void pc_comp_ms_modifiedF_SS(struct tnode_array *tree_array, int idx, int interpolationOrder,
        double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW)
{
    int interpOrderLim = interpolationOrder + 1;
    int pointsPerCluster =  interpOrderLim * interpOrderLim * interpOrderLim;
    int pointsInNode = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClusters = idx * pointsPerCluster;
    int startingIndexInSources = tree_array->ibeg[idx]-1;

    double weights[interpOrderLim], dj[interpOrderLim];
    double nodeX[interpOrderLim], nodeY[interpOrderLim], nodeZ[interpOrderLim];
    
    double *modifiedF, *modifiedF2;
    make_vector(modifiedF, pointsInNode);
    make_vector(modifiedF2, pointsInNode);

    int *exactIndX, *exactIndY, *exactIndZ;
    make_vector(exactIndX, pointsInNode);
    make_vector(exactIndY, pointsInNode);
    make_vector(exactIndZ, pointsInNode);


    double x0 = tree_array->x_min[idx];  // 1e-15 fails for large meshes, mysteriously.
    double x1 = tree_array->x_max[idx];
    double y0 = tree_array->y_min[idx];
    double y1 = tree_array->y_max[idx];
    double z0 = tree_array->z_min[idx];
    double z1 = tree_array->z_max[idx];


#ifdef OPENACC_ENABLED
    int streamID = rand() % 3;
    #pragma acc kernels async(streamID) present(tt, xS, yS, zS, qS, wS, \
                                                clusterX, clusterY, clusterZ, clusterQ, clusterW) \
                       create(modifiedF[0:pointsInNode], modifiedF2[0:pointsInNode], exactIndX[0:pointsInNode], \
                              exactIndY[0:pointsInNode], exactIndZ[0:pointsInNode], \
                              nodeX[0:interpOrderLim], nodeY[0:interpOrderLim], nodeZ[0:interpOrderLim], \
                              weights[0:interpOrderLim], dj[0:interpOrderLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < pointsInNode; j++) {
        modifiedF[j] = qS[startingIndexInSources + j] * wS[startingIndexInSources + j];
        modifiedF2[j] = wS[startingIndexInSources + j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interpOrderLim; i++) {
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interpOrderLim; j++) {
        dj[j] = 1.0;
        if (j == 0) dj[j] = 0.5;
        if (j == interpolationOrder) dj[j] = 0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interpOrderLim; j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

    // Compute modified f values
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < pointsInNode; i++) { // loop through the source points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double sx = xS[startingIndexInSources+i];
        double sy = yS[startingIndexInSources+i];
        double sz = zS[startingIndexInSources+i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:sumX) reduction(+:sumY) reduction(+:sumZ)
#endif
        for (int j = 0; j < interpOrderLim; j++) {  // loop through the degree

            double cx = sx - nodeX[j];
            double cy = sy - nodeY[j];
            double cz = sz - nodeZ[j];

            if (fabs(cx) < DBL_MIN) exactIndX[i] = j;
            if (fabs(cy) < DBL_MIN) exactIndY[i] = j;
            if (fabs(cz) < DBL_MIN) exactIndZ[i] = j;

            // Increment the sums
            double w = weights[j];
            sumX += w / (cx);
            sumY += w / (cy);
            sumZ += w / (cz);

        }

        double denominator = 1.0;
        if (exactIndX[i] == -1) denominator *= sumX;
        if (exactIndY[i] == -1) denominator *= sumY;
        if (exactIndZ[i] == -1) denominator *= sumZ;

        modifiedF[i] /= denominator;
        modifiedF2[i] /= denominator;

    }


#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < pointsPerCluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point
        // compute k1, k2, k3 from j
        int k1 = j % interpOrderLim;
        int kk = (j-k1) / interpOrderLim;
        int k2 = kk % interpOrderLim;
        kk = kk - k2;
        int k3 = kk / interpOrderLim;

        double cz = nodeZ[k3];
        double w3 = weights[k3];

        double cy = nodeY[k2];
        double w2 = weights[k2];

        double cx = nodeX[k1];
        double w1 = weights[k1];

        // Fill cluster X, Y, and Z arrays
        clusterX[startingIndexInClusters + j] = cx;
        clusterY[startingIndexInClusters + j] = cy;
        clusterZ[startingIndexInClusters + j] = cz;

        // Increment cluster Q array
        double temp = 0.0, temp2 = 0.0;

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp) reduction(+:temp2)
#endif
        for (int i = 0; i < pointsInNode; i++) {  // loop over source points
            double sx = xS[startingIndexInSources + i];
            double sy = yS[startingIndexInSources + i];
            double sz = zS[startingIndexInSources + i];

            double numerator = 1.0;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
            if (exactIndX[i] == -1) {
                numerator *=  w1 / (sx - cx);
            } else {
                if (exactIndX[i] != k1) numerator *= 0;
            }

            if (exactIndY[i] == -1) {
                numerator *=  w2 / (sy - cy);
            } else {
                if (exactIndY[i] != k2) numerator *= 0;
            }

            if (exactIndZ[i] == -1) {
                numerator *=  w3 / (sz - cz);
            } else {
                if (exactIndZ[i] != k3) numerator *= 0;
            }

            temp += numerator * modifiedF[i];
            temp2 += numerator * modifiedF2[i];
        }
        
        clusterQ[startingIndexInClusters + j] += temp;
        clusterW[startingIndexInClusters + j] += temp2;
    }
    
#ifdef OPENACC_ENABLED
    } // end acc kernels region
#endif

    free_vector(modifiedF);
    free_vector(modifiedF2);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}



void pc_comp_ms_modifiedF_hermite(struct tnode_array *tree_array, int idx, int interpolationOrder,
        int totalNumberInterpolationPoints, double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW)
{

    int interpOrderLim = interpolationOrder + 1;
    int interpolationPointsPerCluster =  interpOrderLim * interpOrderLim * interpOrderLim;
    int sourcePointsInCluster = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInSourcesArray = tree_array->ibeg[idx] - 1;


    double dj[interpOrderLim], wx[interpOrderLim], wy[interpOrderLim], wz[interpOrderLim];
    double nodeX[interpOrderLim], nodeY[interpOrderLim], nodeZ[interpOrderLim];

    double *modifiedF;
    make_vector(modifiedF, sourcePointsInCluster);

    int *exactIndX, *exactIndY, *exactIndZ;
    make_vector(exactIndX, sourcePointsInCluster);
    make_vector(exactIndY, sourcePointsInCluster);
    make_vector(exactIndZ, sourcePointsInCluster);


    // Set the bounding box.
    double x0 = tree_array->x_min[idx];
    double x1 = tree_array->x_max[idx];
    double y0 = tree_array->y_min[idx];
    double y1 = tree_array->y_max[idx];
    double z0 = tree_array->z_min[idx];
    double z1 = tree_array->z_max[idx];


#ifdef OPENACC_ENABLED
    int streamID = rand() % 3;
    #pragma acc kernels async(streamID) present(tt, ww, xS, yS, zS, qS, wS, \
                                                clusterX, clusterY, clusterZ, clusterQ) \
                       create(modifiedF[0:sourcePointsInCluster], exactIndX[0:sourcePointsInCluster], \
                              exactIndY[0:sourcePointsInCluster], exactIndZ[0:sourcePointsInCluster], \
                              nodeX[0:interpOrderLim], nodeY[0:interpOrderLim], nodeZ[0:interpOrderLim], \
                              dj[0:interpOrderLim], wx[0:interpOrderLim], \
                              wy[0:interpOrderLim], wz[0:interpOrderLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < sourcePointsInCluster; j++) {
        modifiedF[j] = qS[startingIndexInSourcesArray + j] * wS[startingIndexInSourcesArray + j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interpOrderLim; i++) {
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interpOrderLim; j++) {
        dj[j] = 1.0;
        wx[j] = -4.0 * ww[j] / (x1 - x0);
        wy[j] = -4.0 * ww[j] / (y1 - y0);
        wz[j] = -4.0 * ww[j] / (z1 - z0);
    }
    dj[0] = 0.25;
    dj[interpolationOrder] = 0.25;



#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < sourcePointsInCluster; i++) { // loop through the source points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double sx = xS[startingIndexInSourcesArray + i];
        double sy = yS[startingIndexInSourcesArray + i];
        double sz = zS[startingIndexInSourcesArray + i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:sumX) reduction(+:sumY) reduction(+:sumZ)
#endif
        for (int j = 0; j < interpOrderLim; j++) {  // loop through the degree

            double dx = sx - nodeX[j];
            double dy = sy - nodeY[j];
            double dz = sz - nodeZ[j];

            if (fabs(dx) < DBL_MIN) exactIndX[i] = j;
            if (fabs(dy) < DBL_MIN) exactIndY[i] = j;
            if (fabs(dz) < DBL_MIN) exactIndZ[i] = j;

            // Increment the sums
            sumX += dj[j] / (dx*dx) + wx[j] / dx;
            sumY += dj[j] / (dy*dy) + wy[j] / dy;
            sumZ += dj[j] / (dz*dz) + wz[j] / dz;

        }

        double denominator = 1.0;
        if (exactIndX[i] == -1) denominator *= sumX;
        if (exactIndY[i] == -1) denominator *= sumY;
        if (exactIndZ[i] == -1) denominator *= sumZ;

        modifiedF[i] /= denominator;
    }


#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interpolationPointsPerCluster; j++) {
        // compute k1, k2, k3 from j
        int k1 = j % interpOrderLim;
        int kk = (j-k1) / interpOrderLim;
        int k2 = kk % interpOrderLim;
        kk = kk - k2;
        int k3 = kk / interpOrderLim;

        double cz = nodeZ[k3];
        double cy = nodeY[k2];
        double cx = nodeX[k1];

        // Fill cluster X, Y, and Z arrays
        int interpolationPointIndex = startingIndexInClustersArray + j;
        clusterX[interpolationPointIndex] = cx;
        clusterY[interpolationPointIndex] = cy;
        clusterZ[interpolationPointIndex] = cz;


        // Increment cluster Q array
        double temp0 = 0.0, temp1 = 0.0, temp2 = 0.0, temp3 = 0.0;
        double temp4 = 0.0, temp5 = 0.0, temp6 = 0.0, temp7 = 0.0;
         
#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp0) reduction(+:temp1) reduction(+:temp2) \
                                     reduction(+:temp3) reduction(+:temp4) reduction(+:temp5) \
                                     reduction(+:temp6) reduction(+:temp7)
#endif
        for (int i = 0; i < sourcePointsInCluster; i++) {  // loop over source points
        
            int sourcePointIndex = startingIndexInSourcesArray + i;
            double dx = xS[sourcePointIndex] - cx;
            double dy = yS[sourcePointIndex] - cy;
            double dz = zS[sourcePointIndex] - cz;

            double numerator0 = 1.0, numerator1 = 1.0, numerator2 = 1.0, numerator3 = 1.0;
            double numerator4 = 1.0, numerator5 = 1.0, numerator6 = 1.0, numerator7 = 1.0;

            double Ax = dj[k1] / (dx*dx) + wx[k1] / dx;
            double Ay = dj[k2] / (dy*dy) + wy[k2] / dy;
            double Az = dj[k3] / (dz*dz) + wz[k3] / dz;
            double Bx = dj[k1] / dx;
            double By = dj[k2] / dy;
            double Bz = dj[k3] / dz;


            if (exactIndX[i] == -1) {
                numerator0 *=  Ax;                     // Aaa

                numerator1 *=  Bx;                     // Baa
                numerator2 *=  Ax;                     // Aba
                numerator3 *=  Ax;                     // Aab

                numerator4 *=  Bx;                     // Bba
                numerator5 *=  Ax;                     // Abb
                numerator6 *=  Bx;                     // Bab

                numerator7 *=  Bx;                     // Bbb

            } else {
                if (exactIndX[i] != k1) {
                    numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0;
                    numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                } else {
                    numerator1 *= 0; numerator4 *= 0; numerator6 *= 0; numerator7 *= 0;
                }
            }

            if (exactIndY[i] == -1) {
                numerator0 *=  Ay;                    // aAa

                numerator1 *=  Ay;                     // bAa
                numerator2 *=  By;                     // aBa
                numerator3 *=  Ay;                     // aAb

                numerator4 *=  By;                     // bBa
                numerator5 *=  By;                     // aBb
                numerator6 *=  Ay;                     // bAb

                numerator7 *=  By;                     // bBb

            } else {
                if (exactIndY[i] != k2) {
                    numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0;
                    numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                }  else {
                    numerator2 *= 0; numerator4 *= 0; numerator5 *= 0; numerator7 *= 0;
                }

            }

            if (exactIndZ[i] == -1) {
                numerator0 *=  Az;                    // aaA

                numerator1 *=  Az;                    // baA
                numerator2 *=  Az;                    // abA
                numerator3 *=  Bz;                    // aaB

                numerator4 *=  Az;                    // bbA
                numerator5 *=  Bz;                    // abB
                numerator6 *=  Bz;                    // baB
                
                numerator7 *=  Bz;                    // bbB

            } else {
                if (exactIndZ[i] != k3) {
                    numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0;
                    numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                } else {
                    numerator3 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                }
            }

            temp0 += numerator0 * modifiedF[i];
            temp1 += numerator1 * modifiedF[i];
            temp2 += numerator2 * modifiedF[i];
            temp3 += numerator3 * modifiedF[i];
            temp4 += numerator4 * modifiedF[i];
            temp5 += numerator5 * modifiedF[i];
            temp6 += numerator6 * modifiedF[i];
            temp7 += numerator7 * modifiedF[i];

        }

        clusterQ[0 * totalNumberInterpolationPoints + interpolationPointIndex] += temp0;
        clusterQ[1 * totalNumberInterpolationPoints + interpolationPointIndex] += temp1;
        clusterQ[2 * totalNumberInterpolationPoints + interpolationPointIndex] += temp2;
        clusterQ[3 * totalNumberInterpolationPoints + interpolationPointIndex] += temp3;
        clusterQ[4 * totalNumberInterpolationPoints + interpolationPointIndex] += temp4;
        clusterQ[5 * totalNumberInterpolationPoints + interpolationPointIndex] += temp5;
        clusterQ[6 * totalNumberInterpolationPoints + interpolationPointIndex] += temp6;
        clusterQ[7 * totalNumberInterpolationPoints + interpolationPointIndex] += temp7;

    }

#ifdef OPENACC_ENABLED
    } // end acc kernels region
#endif

    free_vector(modifiedF);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}



void pc_comp_ms_modifiedF_hermite_SS(struct tnode_array *tree_array, int idx, int interpolationOrder,
        int totalNumberInterpolationPoints, double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW)
{

    int interpOrderLim = interpolationOrder + 1;
    int interpolationPointsPerCluster =  interpOrderLim * interpOrderLim * interpOrderLim;
    int sourcePointsInCluster = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInSourcesArray = tree_array->ibeg[idx] - 1;


    double dj[interpOrderLim], wx[interpOrderLim], wy[interpOrderLim], wz[interpOrderLim];
    double nodeX[interpOrderLim], nodeY[interpOrderLim], nodeZ[interpOrderLim];

    double *modifiedF, *modifiedF2;
    make_vector(modifiedF, sourcePointsInCluster);
    make_vector(modifiedF2, sourcePointsInCluster);

    int *exactIndX, *exactIndY, *exactIndZ;
    make_vector(exactIndX, sourcePointsInCluster);
    make_vector(exactIndY, sourcePointsInCluster);
    make_vector(exactIndZ, sourcePointsInCluster);


    // Set the bounding box.
    double x0 = tree_array->x_min[idx];
    double x1 = tree_array->x_max[idx];
    double y0 = tree_array->y_min[idx];
    double y1 = tree_array->y_max[idx];
    double z0 = tree_array->z_min[idx];
    double z1 = tree_array->z_max[idx];


#ifdef OPENACC_ENABLED
    int streamID = rand() % 3;
    #pragma acc kernels async(streamID) present(tt, ww, xS, yS, zS, qS, wS, \
                                                clusterX, clusterY, clusterZ, clusterQ) \
                       create(modifiedF[0:sourcePointsInCluster], modifiedF2[0:sourcePointsInCluster], \
                              exactIndX[0:sourcePointsInCluster], exactIndY[0:sourcePointsInCluster], \
                              exactIndZ[0:sourcePointsInCluster], \
                              nodeX[0:interpOrderLim], nodeY[0:interpOrderLim], nodeZ[0:interpOrderLim], \
                              dj[0:interpOrderLim], wx[0:interpOrderLim], \
                              wy[0:interpOrderLim], wz[0:interpOrderLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < sourcePointsInCluster; j++) {
        modifiedF[j] = qS[startingIndexInSourcesArray + j] * wS[startingIndexInSourcesArray + j];
        modifiedF2[j] = wS[startingIndexInSourcesArray + j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interpOrderLim; i++) {
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interpOrderLim; j++) {
        dj[j] = 1.0;
        wx[j] = -4.0 * ww[j] / (x1 - x0);
        wy[j] = -4.0 * ww[j] / (y1 - y0);
        wz[j] = -4.0 * ww[j] / (z1 - z0);
    }
    dj[0] = 0.25;
    dj[interpolationOrder] = 0.25;



#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < sourcePointsInCluster; i++) { // loop through the source points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double sx = xS[startingIndexInSourcesArray + i];
        double sy = yS[startingIndexInSourcesArray + i];
        double sz = zS[startingIndexInSourcesArray + i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:sumX) reduction(+:sumY) reduction(+:sumZ)
#endif
        for (int j = 0; j < interpOrderLim; j++) {  // loop through the degree

            double dx = sx - nodeX[j];
            double dy = sy - nodeY[j];
            double dz = sz - nodeZ[j];

            if (fabs(dx) < DBL_MIN) exactIndX[i] = j;
            if (fabs(dy) < DBL_MIN) exactIndY[i] = j;
            if (fabs(dz) < DBL_MIN) exactIndZ[i] = j;

            // Increment the sums
            sumX += dj[j] / (dx*dx) + wx[j] / dx;
            sumY += dj[j] / (dy*dy) + wy[j] / dy;
            sumZ += dj[j] / (dz*dz) + wz[j] / dz;

        }

        double denominator = 1.0;
        if (exactIndX[i] == -1) denominator *= sumX;
        if (exactIndY[i] == -1) denominator *= sumY;
        if (exactIndZ[i] == -1) denominator *= sumZ;

        modifiedF[i] /= denominator;
        modifiedF2[i] /= denominator;
    } 


#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interpolationPointsPerCluster; j++) {
        // compute k1, k2, k3 from j
        int k1 = j % interpOrderLim;
        int kk = (j-k1) / interpOrderLim;
        int k2 = kk % interpOrderLim;
        kk = kk - k2;
        int k3 = kk / interpOrderLim;

        double cz = nodeZ[k3];
        double cy = nodeY[k2];
        double cx = nodeX[k1];

        // Fill cluster X, Y, and Z arrays
        int interpolationPointIndex = startingIndexInClustersArray + j;
        clusterX[interpolationPointIndex] = cx;
        clusterY[interpolationPointIndex] = cy;
        clusterZ[interpolationPointIndex] = cz;


        // Increment cluster Q array
        double tempq0 = 0.0, tempq1 = 0.0, tempq2 = 0.0, tempq3 = 0.0;
        double tempq4 = 0.0, tempq5 = 0.0, tempq6 = 0.0, tempq7 = 0.0;

        double tempw0 = 0.0, tempw1 = 0.0, tempw2 = 0.0, tempw3 = 0.0;
        double tempw4 = 0.0, tempw5 = 0.0, tempw6 = 0.0, tempw7 = 0.0;
         
#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:tempq0) reduction(+:tempq1) reduction(+:tempq2) \
                                     reduction(+:tempq3) reduction(+:tempq4) reduction(+:tempq5) \
                                     reduction(+:tempq6) reduction(+:tempq7) \
                                     reduction(+:tempw0) reduction(+:tempw1) reduction(+:tempw2) \
                                     reduction(+:tempw3) reduction(+:tempw4) reduction(+:tempw5) \
                                     reduction(+:tempw6) reduction(+:tempw7)
#endif
        for (int i = 0; i < sourcePointsInCluster; i++) {  // loop over source points
         
             int sourcePointIndex = startingIndexInSourcesArray + i;
             double dx = xS[sourcePointIndex] - cx;
             double dy = yS[sourcePointIndex] - cy;
             double dz = zS[sourcePointIndex] - cz;

             double numerator0 = 1.0, numerator1 = 1.0, numerator2 = 1.0, numerator3 = 1.0;
             double numerator4 = 1.0, numerator5 = 1.0, numerator6 = 1.0, numerator7 = 1.0;

             double Ax = dj[k1] / (dx*dx) + wx[k1] / dx;
             double Ay = dj[k2] / (dy*dy) + wy[k2] / dy;
             double Az = dj[k3] / (dz*dz) + wz[k3] / dz;
             double Bx = dj[k1] / dx;
             double By = dj[k2] / dy;
             double Bz = dj[k3] / dz;


             if (exactIndX[i] == -1) {
                 numerator0 *=  Ax;                     // Aaa

                 numerator1 *=  Bx;                     // Baa
                 numerator2 *=  Ax;                     // Aba
                 numerator3 *=  Ax;                     // Aab

                 numerator4 *=  Bx;                     // Bba
                 numerator5 *=  Ax;                     // Abb
                 numerator6 *=  Bx;                     // Bab

                 numerator7 *=  Bx;                     // Bbb

             } else {
                 if (exactIndX[i] != k1) {
                     numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0;
                     numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                 } else {
                     numerator1 *= 0; numerator4 *= 0; numerator6 *= 0; numerator7 *= 0;
                 }
             }

             if (exactIndY[i] == -1) {
                 numerator0 *=  Ay;                    // aAa

                 numerator1 *=  Ay;                     // bAa
                 numerator2 *=  By;                     // aBa
                 numerator3 *=  Ay;                     // aAb

                 numerator4 *=  By;                     // bBa
                 numerator5 *=  By;                     // aBb
                 numerator6 *=  Ay;                     // bAb

                 numerator7 *=  By;                     // bBb

             } else {
                 if (exactIndY[i] != k2) {
                     numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0;
                     numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                 }  else {
                     numerator2 *= 0; numerator4 *= 0; numerator5 *= 0; numerator7 *= 0;
                 }

             }

             if (exactIndZ[i] == -1) {
                 numerator0 *=  Az;                    // aaA

                 numerator1 *=  Az;                    // baA
                 numerator2 *=  Az;                    // abA
                 numerator3 *=  Bz;                    // aaB

                 numerator4 *=  Az;                    // bbA
                 numerator5 *=  Bz;                    // abB
                 numerator6 *=  Bz;                    // baB
                 
                 numerator7 *=  Bz;                    // bbB

             } else {
                 if (exactIndZ[i] != k3) {
                     numerator0 *= 0; numerator1 *= 0; numerator2 *= 0; numerator3 *= 0;
                     numerator4 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                 } else {
                     numerator3 *= 0; numerator5 *= 0; numerator6 *= 0; numerator7 *= 0;
                 }
             }

             tempq0 += numerator0 * modifiedF[i];
             tempq1 += numerator1 * modifiedF[i];
             tempq2 += numerator2 * modifiedF[i];
             tempq3 += numerator3 * modifiedF[i];
             tempq4 += numerator4 * modifiedF[i];
             tempq5 += numerator5 * modifiedF[i];
             tempq6 += numerator6 * modifiedF[i];
             tempq7 += numerator7 * modifiedF[i];

             tempw0 += numerator0 * modifiedF2[i];
             tempw1 += numerator1 * modifiedF2[i];
             tempw2 += numerator2 * modifiedF2[i];
             tempw3 += numerator3 * modifiedF2[i];
             tempw4 += numerator4 * modifiedF2[i];
             tempw5 += numerator5 * modifiedF2[i];
             tempw6 += numerator6 * modifiedF2[i];
             tempw7 += numerator7 * modifiedF2[i];

         }

         clusterQ[0 * totalNumberInterpolationPoints + interpolationPointIndex] += tempq0;
         clusterQ[1 * totalNumberInterpolationPoints + interpolationPointIndex] += tempq1;
         clusterQ[2 * totalNumberInterpolationPoints + interpolationPointIndex] += tempq2;
         clusterQ[3 * totalNumberInterpolationPoints + interpolationPointIndex] += tempq3;
         clusterQ[4 * totalNumberInterpolationPoints + interpolationPointIndex] += tempq4;
         clusterQ[5 * totalNumberInterpolationPoints + interpolationPointIndex] += tempq5;
         clusterQ[6 * totalNumberInterpolationPoints + interpolationPointIndex] += tempq6;
         clusterQ[7 * totalNumberInterpolationPoints + interpolationPointIndex] += tempq7;

         clusterW[0 * totalNumberInterpolationPoints + interpolationPointIndex] += tempw0;
         clusterW[1 * totalNumberInterpolationPoints + interpolationPointIndex] += tempw1;
         clusterW[2 * totalNumberInterpolationPoints + interpolationPointIndex] += tempw2;
         clusterW[3 * totalNumberInterpolationPoints + interpolationPointIndex] += tempw3;
         clusterW[4 * totalNumberInterpolationPoints + interpolationPointIndex] += tempw4;
         clusterW[5 * totalNumberInterpolationPoints + interpolationPointIndex] += tempw5;
         clusterW[6 * totalNumberInterpolationPoints + interpolationPointIndex] += tempw6;
         clusterW[7 * totalNumberInterpolationPoints + interpolationPointIndex] += tempw7;

    }

#ifdef OPENACC_ENABLED
    } // end acc kernels region
#endif

    free_vector(modifiedF);
    free_vector(modifiedF2);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}
