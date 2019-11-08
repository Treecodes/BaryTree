
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

void pc_comp_ms_modifiedF(struct tnode_array * tree_array, int idx, int interpolationOrder,
                          double *xS, double *yS, double *zS, double *qS, double *wS,
                          double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW);

void pc_comp_ms_modifiedF_SS(struct tnode_array * tree_array, int idx, int interpolationOrder,
                             double *xS, double *yS, double *zS, double *qS, double *wS,
                             double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW);


void pc_comp_ms_modifiedF_hermite(struct tnode *p, double *xS, double *yS, double *zS, double *qS, double *wS,
                                  double *clusterX, double *clusterY, double *clusterZ, double *clusterQ,
                                  double *clusterMx, double *clusterMy, double *clusterMz,
                                  double *clusterMxy, double *clusterMyz, double *clusterMzx, double *clusterMxyz);

void pc_comp_ms_modifiedF_hermite_SS(struct tnode *p, double *xS, double *yS, double *zS, double *qS, double *wS,
                                     double *clusterX, double *clusterY, double *clusterZ, double *clusterQ,
                                     double *clusterMx, double *clusterMy, double *clusterMz,
                                     double *clusterMxy, double *clusterMyz, double *clusterMzx, double *clusterMxyz,
                                     double *clusterW, double *clusterWx, double *clusterWy, double *clusterWz,
                                     double *clusterWxy, double *clusterWyz, double *clusterWzx, double *clusterWxyz);



void Clusters_PC_SetupLagrange(struct particles *clusters, struct particles *sources, struct tnode *troot,
                               int interpolationOrder, struct tnode_array *tree_array, char *singularityHandling)
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
    }

    if (strcmp(singularityHandling, "skipping") == 0) {
        for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->w[i] = 1.0;

    } else if (strcmp(singularityHandling, "subtraction") == 0) {
        for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->w[i] = 0.0;
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

    if (strcmp(singularityHandling, "skipping") == 0) {
        for (int i = 0; i < tree_numnodes; i++)
            pc_comp_ms_modifiedF(tree_array, i, interpolationOrder, xS, yS, zS, qS, wS, xC, yC, zC, qC, wC);

    } else if (strcmp(singularityHandling, "subtraction") == 0) {
        for (int i = 0; i < tree_numnodes; i++)
            pc_comp_ms_modifiedF_SS(tree_array, i, interpolationOrder, xS, yS, zS, qS, wS, xC, yC, zC, qC, wC);
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


void pc_comp_ms_modifiedF_SS(struct tnode_array *tree_array, int idx, int interpolationOrder,
        double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW)
{
    int i,j,k;

    int interpOrderLim = interpolationOrder + 1;
    int pointsPerCluster =  interpOrderLim * interpOrderLim * interpOrderLim;
    int pointsInNode = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClusters = idx * pointsPerCluster;
    int startingIndexInSources = tree_array->ibeg[idx]-1;

    double x0, x1, y0, y1, z0, z1;  // bounding box

    double weights[interpOrderLim], dj[interpOrderLim];
    double *modifiedF, *modifiedF2;
    make_vector(modifiedF,pointsInNode);
    make_vector(modifiedF2,pointsInNode);

    double nodeX[interpOrderLim], nodeY[interpOrderLim], nodeZ[interpOrderLim];

    int *exactIndX, *exactIndY, *exactIndZ;
    make_vector(exactIndX, pointsInNode);
    make_vector(exactIndY, pointsInNode);
    make_vector(exactIndZ, pointsInNode);

    x0 = tree_array->x_min[idx];  // 1e-15 fails for large meshes, mysteriously.
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
    #pragma acc kernels async(streamID) present(xS, yS, zS, qS, wS, clusterX, clusterY, clusterZ, clusterQ, clusterW,tt) \
    create(modifiedF[0:pointsInNode],modifiedF2[0:pointsInNode],exactIndX[0:pointsInNode],exactIndY[0:pointsInNode],exactIndZ[0:pointsInNode], \
            nodeX[0:interpOrderLim],nodeY[0:interpOrderLim],nodeZ[0:interpOrderLim], \
            weights[0:interpOrderLim],dj[0:interpOrderLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < pointsInNode; j++) {
        modifiedF[j] = qS[startingIndexInSources+j] * wS[startingIndexInSources+j];
//        modifiedF[j] = qS[startingIndexInSources+j];
        modifiedF2[j] = wS[startingIndexInSources+j];
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

        sumX = 0.0;
        sumY = 0.0;
        sumZ = 0.0;

        sx = xS[startingIndexInSources+i];
        sy = yS[startingIndexInSources+i];
        sz = zS[startingIndexInSources+i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (int j = 0; j < interpOrderLim; j++) {  // loop through the degree

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
        modifiedF2[i] /= denominator;

    }

    // Compute moments for each interpolation point
    double numerator, xn, yn, zn, temp, temp2;
    int k1, k2, k3, kk;
    double w1,w2,w3;

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < pointsPerCluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point
        // compute k1, k2, k3 from j
        k1 = j%interpOrderLim;
        kk = (j-k1)/interpOrderLim;
        k2 = kk%interpOrderLim;
        kk = kk - k2;
        k3 = kk / interpOrderLim;

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
        temp2=0.0;
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
            temp2 += numerator*modifiedF2[i];
        }
        clusterQ[startingIndexInClusters + j] += temp;
		clusterW[startingIndexInClusters + j] += temp2;
    }
#ifdef OPENACC_ENABLED
    }
#endif

    free_vector(modifiedF);
    free_vector(modifiedF2);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}
