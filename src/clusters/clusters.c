#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../utilities/enums.h"

#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"

#include "struct_clusters.h"
#include "clusters.h"


static void pc_comp_ms_modifiedF(const struct Tree *tree, int idx, int interpolationOrder,
                double *xS, double *yS, double *zS, double *qS,
                double *clusterX, double *clusterY, double *clusterZ, double *clusterQ);

static void cp_comp_interp(const struct Tree *tree, int idx, int interpolationOrder,
                double *clusterX, double *clusterY, double *clusterZ);



void Clusters_Sources_Construct(struct Clusters **clusters_addr, const struct Particles *sources,
                                const struct Tree *tree, const struct RunParams *run_params)
{
    *clusters_addr = malloc(sizeof(struct Clusters));
    struct Clusters *clusters = *clusters_addr;
    
    APPROXIMATION approximation = run_params->approximation;

    int tree_numnodes = tree->numnodes;
    int totalNumberSourcePoints = sources->num;
    
    int interpolationOrder = run_params->interp_order;
    int interpolationPointsPerCluster = run_params->interp_pts_per_cluster;

    int totalNumberInterpolationPoints  = tree_numnodes * interpolationPointsPerCluster;
    int totalNumberInterpolationCharges = tree_numnodes * run_params->interp_charges_per_cluster;

    clusters->x = NULL;
    clusters->y = NULL;
    clusters->z = NULL;
    clusters->q = NULL;
    
    make_vector(clusters->x, totalNumberInterpolationPoints);
    make_vector(clusters->y, totalNumberInterpolationPoints);
    make_vector(clusters->z, totalNumberInterpolationPoints);
    make_vector(clusters->q, totalNumberInterpolationCharges);

    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->x[i] = 0.0;
    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->y[i] = 0.0;
    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->z[i] = 0.0;
    for (int i = 0; i < totalNumberInterpolationCharges; i++) clusters->q[i] = 0.0;
        
    clusters->num = totalNumberInterpolationPoints;
    clusters->num_charges = totalNumberInterpolationCharges;

    double *xS = sources->x;
    double *yS = sources->y;
    double *zS = sources->z;
    double *qS = sources->q;

    double *xC = clusters->x;
    double *yC = clusters->y;
    double *zC = clusters->z;
    double *qC = clusters->q;


#ifdef OPENACC_ENABLED
    #pragma acc data copyin(xS[0:totalNumberSourcePoints], yS[0:totalNumberSourcePoints], \
                            zS[0:totalNumberSourcePoints], qS[0:totalNumberSourcePoints]) \
                       copy(xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], \
                            zC[0:totalNumberInterpolationPoints], qC[0:totalNumberInterpolationCharges])
    {
#endif

    if (approximation == LAGRANGE) {
        for (int i = 0; i < tree_numnodes; i++)
            pc_comp_ms_modifiedF(tree, i, interpolationOrder, xS, yS, zS, qS, xC, yC, zC, qC);

    } else {
        exit(1);
    }
    
#ifdef OPENACC_ENABLED
    #pragma acc wait
    } // end ACC DATA REGION
#endif

    return;
}




void Clusters_Targets_Construct(struct Clusters **clusters_addr, const struct Tree *tree,
                                const struct RunParams *run_params)
{
    *clusters_addr = malloc(sizeof(struct Clusters));
    struct Clusters *clusters = *clusters_addr;

    int tree_numnodes = tree->numnodes;

    int interpolationOrder = run_params->interp_order;
    int interpOrderLim = interpolationOrder + 1;

    int totalNumberInterpolationPoints  = tree_numnodes * interpOrderLim;
    int totalNumberInterpolationCharges = tree_numnodes * run_params->interp_charges_per_cluster;

    clusters->x = NULL;
    clusters->y = NULL;
    clusters->z = NULL;
    clusters->q = NULL;

    make_vector(clusters->x, totalNumberInterpolationPoints);
    make_vector(clusters->y, totalNumberInterpolationPoints);
    make_vector(clusters->z, totalNumberInterpolationPoints);
    make_vector(clusters->q, totalNumberInterpolationCharges);

    clusters->num         = totalNumberInterpolationPoints;
    clusters->num_charges = totalNumberInterpolationCharges;

    double *xC = clusters->x;
    double *yC = clusters->y;
    double *zC = clusters->z;
    double *qC = clusters->q;


#ifdef OPENACC_ENABLED
    #pragma acc enter data create(xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], \
                             zC[0:totalNumberInterpolationPoints], qC[0:totalNumberInterpolationCharges])
    {
#endif

    for (int i = 0; i < tree_numnodes; i++) {
        if ((tree->used_children[i] > 0 && (tree->used[i] == 1 || tree->used_parent[i] == 1))
          || tree->used_leaf[i] == 1) {
            cp_comp_interp(tree, i, interpolationOrder, xC, yC, zC);
        }
    }
    
#ifdef OPENACC_ENABLED
    #pragma acc wait
    } // end ACC DATA REGION
#endif

    return;
}




void Clusters_Alloc(struct Clusters **clusters_addr, int length, const struct RunParams *run_params)
{
    *clusters_addr = malloc(sizeof(struct Clusters));
    struct Clusters *clusters = *clusters_addr;

    APPROXIMATION approximation = run_params->approximation;
    
    clusters->num = length;
    clusters->num_charges = length;
    
    clusters->x = NULL;
    clusters->y = NULL;
    clusters->z = NULL;
    clusters->q = NULL;

    if (approximation == HERMITE)
        clusters->num_charges *= 8;

    if (clusters->num > 0) {
        make_vector(clusters->x, clusters->num);
        make_vector(clusters->y, clusters->num);
        make_vector(clusters->z, clusters->num);
        make_vector(clusters->q, clusters->num_charges);
    }

    return;
}   /* END of function allocate_cluster */




void Clusters_Free(struct Clusters **clusters_addr)
{
    struct Clusters *clusters = *clusters_addr;
    
    if (clusters != NULL) {
#ifdef OPENACC_ENABLED
        #pragma acc exit data delete(clusters->x, clusters->y, clusters->z, clusters->q)
#endif
        if (clusters->x != NULL) {
            free_vector(clusters->x);
        }
        if (clusters->y != NULL) {
            free_vector(clusters->y);
        }
        if (clusters->z != NULL) {
            free_vector(clusters->z);
        }
        if (clusters->q != NULL) {
            free_vector(clusters->q);
        }
        free(clusters);
    }
    
    clusters = NULL;

    return;
}   /* END of function Clusters_Free */



/************************************/
/***** LOCAL FUNCTIONS **************/
/************************************/

void pc_comp_ms_modifiedF(const struct Tree *tree, int idx, int interpolationOrder,
        double *xS, double *yS, double *zS, double *qS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ)
{

    int interpOrderLim = interpolationOrder + 1;
    int interpolationPointsPerCluster = interpOrderLim * interpOrderLim * interpOrderLim;
    int sourcePointsInCluster = tree->iend[idx] - tree->ibeg[idx] + 1;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInSourcesArray = tree->ibeg[idx]-1;

    double *weights, *dj, *tt, *nodeX, *nodeY, *nodeZ, *modifiedF;
    int *exactIndX, *exactIndY, *exactIndZ;

    make_vector(weights,   interpOrderLim);
    make_vector(dj,        interpOrderLim);
    make_vector(tt,        interpOrderLim);
    make_vector(nodeX,     interpOrderLim);
    make_vector(nodeY,     interpOrderLim);
    make_vector(nodeZ,     interpOrderLim);
    make_vector(modifiedF, sourcePointsInCluster);
    make_vector(exactIndX, sourcePointsInCluster);
    make_vector(exactIndY, sourcePointsInCluster);
    make_vector(exactIndZ, sourcePointsInCluster);

    double x0 = tree->x_min[idx];
    double x1 = tree->x_max[idx];
    double y0 = tree->y_min[idx];
    double y1 = tree->y_max[idx];
    double z0 = tree->z_min[idx];
    double z1 = tree->z_max[idx];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(xS, yS, zS, qS, clusterX, clusterY, clusterZ, clusterQ) \
                       create(modifiedF[0:sourcePointsInCluster], exactIndX[0:sourcePointsInCluster], \
                              exactIndY[0:sourcePointsInCluster], exactIndZ[0:sourcePointsInCluster], \
                              nodeX[0:interpOrderLim], nodeY[0:interpOrderLim], \
                              nodeZ[0:interpOrderLim], weights[0:interpOrderLim], \
                              dj[0:interpOrderLim], tt[0:interpOrderLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < sourcePointsInCluster; j++) {
        modifiedF[j] = qS[startingIndexInSourcesArray + j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interpOrderLim; i++) {
        tt[i] = cos(i * M_PI / interpolationOrder);
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

    free_vector(weights);
    free_vector(dj);
    free_vector(tt);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);
    free_vector(modifiedF);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}



void cp_comp_interp(const struct Tree *tree, int idx, int interpolationOrder,
                    double *clusterX, double *clusterY, double *clusterZ)
{

    int interpOrderLim = interpolationOrder + 1;
    int startingIndexInClustersArray = idx * interpOrderLim;

    double x0 = tree->x_min[idx];
    double x1 = tree->x_max[idx];
    double y0 = tree->y_min[idx];
    double y1 = tree->y_max[idx];
    double z0 = tree->z_min[idx];
    double z1 = tree->z_max[idx];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(clusterX, clusterY, clusterZ)
    {
#endif


    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interpOrderLim; i++) {
        double tt = cos(i * M_PI / interpolationOrder);
        clusterX[startingIndexInClustersArray + i] = x0 + (tt + 1.0)/2.0 * (x1 - x0);
        clusterY[startingIndexInClustersArray + i] = y0 + (tt + 1.0)/2.0 * (y1 - y0);
        clusterZ[startingIndexInClustersArray + i] = z0 + (tt + 1.0)/2.0 * (z1 - z0);
    }

#ifdef OPENACC_ENABLED
    } //end acc kernels region
#endif

    return;
}
