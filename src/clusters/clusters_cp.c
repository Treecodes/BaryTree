#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../utilities/enums.h"

#include "../globvars.h"

#include "../tree/struct_nodes.h"
#include "../particles/struct_particles.h"

#include "struct_clusters.h"
#include "clusters.h"


static void cp_comp_interp(struct tnode_array *tree_array, int idx, int interpolationOrder,
                           double *clusterX, double *clusterY, double *clusterZ);


void Clusters_CP_Setup(struct clusters **new_clusters,
                       int interpolationOrder, struct tnode_array *tree_array,
                       APPROXIMATION approxName, SINGULARITY singularity)
{
    *new_clusters = malloc(sizeof(struct clusters));
    struct clusters *clusters = *new_clusters;

    int tree_numnodes = tree_array->numnodes;

    int interpOrderLim = interpolationOrder + 1;
    int interpolationPointsPerCluster = interpOrderLim * interpOrderLim * interpOrderLim;

    int totalNumberInterpolationPoints = tree_numnodes * interpolationPointsPerCluster;
    int totalNumberInterpolationCharges = totalNumberInterpolationPoints;
    int totalNumberInterpolationWeights = totalNumberInterpolationPoints; 

    make_vector(clusters->x, totalNumberInterpolationPoints);
    make_vector(clusters->y, totalNumberInterpolationPoints);
    make_vector(clusters->z, totalNumberInterpolationPoints);
    
//    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->x[i] = 0.0;
//    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->y[i] = 0.0;
//    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->z[i] = 0.0;
    
    
    if (approxName == LAGRANGE) {
        make_vector(clusters->q, totalNumberInterpolationCharges);
        make_vector(clusters->w, totalNumberInterpolationWeights);

        for (int i = 0; i < totalNumberInterpolationCharges; i++) clusters->q[i] = 0.0;
        
        if (singularity == SKIPPING) {
            for (int i = 0; i < totalNumberInterpolationWeights; i++) clusters->w[i] = 1.0;

        } else if (singularity == SUBTRACTION) {
            for (int i = 0; i < totalNumberInterpolationWeights; i++) clusters->w[i] = 0.0;

        } else {
            exit(1);
        }
        
    } else if (approxName == HERMITE) {
        totalNumberInterpolationCharges *= 8;

        make_vector(clusters->q, totalNumberInterpolationCharges);
        for (int i = 0; i < totalNumberInterpolationCharges; i++) clusters->q[i] = 0.0;
        
        if (singularity == SKIPPING) {

            make_vector(clusters->w, totalNumberInterpolationWeights);
            for (int i = 0; i < totalNumberInterpolationWeights; i++) clusters->w[i] = 1.0;

        } else if (singularity == SUBTRACTION) {
            totalNumberInterpolationWeights *= 8;

            make_vector(clusters->w, totalNumberInterpolationWeights);
            for (int i = 0; i < totalNumberInterpolationWeights; i++) clusters->w[i] = 0.0;

        } else {
            exit(1);
        }
        
    } else {
        exit(1);
    }

    clusters->num = totalNumberInterpolationPoints;
    clusters->num_charges = totalNumberInterpolationCharges;
    clusters->num_weights = totalNumberInterpolationWeights;

    double *xC = clusters->x;
    double *yC = clusters->y;
    double *zC = clusters->z;

//    printf("First source weight in clusters routine = %f\n", wS[0]);


#ifdef OPENACC_ENABLED
    #pragma acc data copyin(tt[0:interpOrderLim]) \
                     copyout(xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], \
                             zC[0:totalNumberInterpolationPoints])
    {
#endif

    for (int i = 0; i < tree_numnodes; i++)
        cp_comp_interp(tree_array, i, interpolationOrder, xC, yC, zC);
    
#ifdef OPENACC_ENABLED
    #pragma acc wait
    } // end ACC DATA REGION
#endif

    return;
}



/************************************/
/***** LOCAL FUNCTIONS **************/
/************************************/


void cp_comp_interp(struct tnode_array *tree_array, int idx, int interpolationOrder,
                double *clusterX, double *clusterY, double *clusterZ)
{

    int interpOrderLim = interpolationOrder + 1;
    int interpolationPointsPerCluster = interpOrderLim * interpOrderLim * interpOrderLim;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;

    double nodeX[interpOrderLim], nodeY[interpOrderLim], nodeZ[interpOrderLim];

    double x0 = tree_array->x_min[idx];
    double x1 = tree_array->x_max[idx];
    double y0 = tree_array->y_min[idx];
    double y1 = tree_array->y_max[idx];
    double z0 = tree_array->z_min[idx];
    double z1 = tree_array->z_max[idx];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(clusterX, clusterY, clusterZ, tt) \
                       create(nodeX[0:interpOrderLim], nodeY[0:interpOrderLim], \
                              nodeZ[0:interpOrderLim])
    {
#endif


    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interpOrderLim; i++) {
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
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

        // Fill cluster X, Y, and Z arrays
        clusterX[startingIndexInClustersArray + j] = nodeX[k1];
        clusterY[startingIndexInClustersArray + j] = nodeY[k2];
        clusterZ[startingIndexInClustersArray + j] = nodeZ[k3];
    }
#ifdef OPENACC_ENABLED
    } //end acc kernels region
#endif

    return;
}
