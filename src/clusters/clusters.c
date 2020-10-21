#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../utilities/enums.h"

#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"

#include "struct_clusters.h"
#include "clusters.h"


static void pc_comp_ms_modifiedF(const struct Tree *tree, int idx, int interpolationDegree,
                double *xS, double *yS, double *zS, double *qS,
                double *clusterX, double *clusterY, double *clusterZ, double *clusterQ);

static void pc_comp_ms_modifiedF_child_to_parent(const struct Tree *tree, int child_index, int parent_index, int interpolationDegree,
                double *clusterX, double *clusterY, double *clusterZ, double *clusterQ);

static void pc_comp_ms_modifiedF_SS_child_to_parent(const struct Tree *tree, int child_index, int parent_index, int interpolationDegree,
                double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW);

static void pc_comp_ms_modifiedF_SS(const struct Tree *tree, int idx, int interpolationDegree,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW);

static void pc_comp_ms_modifiedF_hermite(const struct Tree *tree, int idx, int interpolationDegree,
                int totalNumberInterpolationPoints,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *clusterX, double *clusterY, double *clusterZ, double *clusterQ);

static void pc_comp_ms_modifiedF_hermite_SS(const struct Tree *tree, int idx, int interpolationDegree,
                int totalNumberInterpolationPoints,
                double *xS, double *yS, double *zS, double *qS, double *wS,
                double *clusterX, double *clusterY, double *clusterZ,
                double *clusterQ, double *clusterW);

static void cp_comp_interp(const struct Tree *tree, int idx, int interpolationDegree,
                double *clusterX, double *clusterY, double *clusterZ);



void Clusters_Sources_Construct(struct Clusters **clusters_addr, const struct Particles *sources,
                                const struct Tree *tree, const struct RunParams *run_params)
{
    *clusters_addr = malloc(sizeof(struct Clusters));
    struct Clusters *clusters = *clusters_addr;
    
    APPROXIMATION approximation = run_params->approximation;
    SINGULARITY singularity = run_params->singularity;

    int tree_numnodes = tree->numnodes;
    int totalNumberSourcePoints = sources->num;
    
    int interpolationDegree = run_params->interp_degree;
    int interpDegreeLim = interpolationDegree + 1;
    int interpolationPointsPerCluster = run_params->interp_pts_per_cluster;

    int totalNumberInterpolationPoints  = tree_numnodes * interpolationPointsPerCluster;
    int totalNumberInterpolationCharges = tree_numnodes * run_params->interp_charges_per_cluster;
    int totalNumberInterpolationWeights = tree_numnodes * run_params->interp_weights_per_cluster;

    clusters->x = NULL;
    clusters->y = NULL;
    clusters->z = NULL;
    clusters->q = NULL;
    clusters->w = NULL;
    
    MPI_Alloc_mem(totalNumberInterpolationPoints*sizeof(double), MPI_INFO_NULL,  &(clusters->x));
    MPI_Alloc_mem(totalNumberInterpolationPoints*sizeof(double), MPI_INFO_NULL,  &(clusters->y));
    MPI_Alloc_mem(totalNumberInterpolationPoints*sizeof(double), MPI_INFO_NULL,  &(clusters->z));
    MPI_Alloc_mem(totalNumberInterpolationCharges*sizeof(double), MPI_INFO_NULL, &(clusters->q));

    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->x[i] = 0.0;
    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->y[i] = 0.0;
    for (int i = 0; i < totalNumberInterpolationPoints; i++) clusters->z[i] = 0.0;
    for (int i = 0; i < totalNumberInterpolationCharges; i++) clusters->q[i] = 0.0;

    if (singularity == SUBTRACTION) {
        MPI_Alloc_mem(totalNumberInterpolationWeights*sizeof(double), MPI_INFO_NULL, &(clusters->w));
        for (int i = 0; i < totalNumberInterpolationWeights; i++) clusters->w[i] = 0.0;
    }

    clusters->num = totalNumberInterpolationPoints;
    clusters->num_charges = totalNumberInterpolationCharges;
    clusters->num_weights = totalNumberInterpolationWeights;

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
    #pragma acc enter data create(xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], \
                                  zC[0:totalNumberInterpolationPoints], qC[0:totalNumberInterpolationCharges])

    #pragma acc kernels present(xC,yC,zC,qC)
    {
        for (int i=0;i<totalNumberInterpolationPoints;i++){
            xC[i]=0.0;
            yC[i]=0.0;
            zC[i]=0.0;
        }

        for (int i=0;i<totalNumberInterpolationCharges;i++){
            qC[i]=0.0;
        }
    }

if (singularity == SUBTRACTION) {
    #pragma acc enter data create(wC[0:totalNumberInterpolationWeights])

    #pragma acc kernels present(wC)
    {
        for (int i=0;i<totalNumberInterpolationWeights;i++){
            wC[i]=0.0;
        }
    }

    }
#endif

    if ((approximation == LAGRANGE) && (singularity == SKIPPING)) {

        // anterpolate from particles to leaf cluster interpolation points
        for (int i = 0; i < tree->leaves_list_num; ++i) {
            int leaf_index = tree->leaves_list[i];
            pc_comp_ms_modifiedF(tree, leaf_index, interpolationDegree, xS, yS, zS, qS, xC, yC, zC, qC);
        }

        // anterpolate up clusters, level by level
        for (int level = tree->max_depth-1; level >= 0; --level) {
            for (int cluster_index = 0; cluster_index < tree->levels_list_num[level]; ++cluster_index) {

                int parent_index = tree->levels_list[level][cluster_index];

                for (int child_counter=0; child_counter<tree->num_children[parent_index]; ++child_counter){

                    int child_index = tree->children[8*parent_index + child_counter];

                    pc_comp_ms_modifiedF_child_to_parent(tree, child_index, parent_index, interpolationDegree, xC, yC, zC, qC);
                }
            }
        }


    } else if ((approximation == LAGRANGE) && (singularity == SUBTRACTION)) {

        // anterpolate from particles to leaf cluster interpolation points
        for (int i = 0; i < tree->leaves_list_num; ++i) {
            int leaf_index = tree->leaves_list[i];
            pc_comp_ms_modifiedF_SS(tree, leaf_index, interpolationDegree, xS, yS, zS, qS, wS, xC, yC, zC, qC, wC);
        }

        // interpolate up clusters, level by level
        for (int level = tree->max_depth-1; level >= 0; --level) {
            for (int cluster_index = 0; cluster_index < tree->levels_list_num[level]; ++cluster_index) {

                int parent_index = tree->levels_list[level][cluster_index];

                for (int child_counter=0; child_counter<tree->num_children[parent_index]; ++child_counter){

                    int child_index = tree->children[8*parent_index + child_counter];

                    pc_comp_ms_modifiedF_SS_child_to_parent(tree, child_index, parent_index, interpolationDegree, xC, yC, zC, qC, wC);
                }
            }
        }

    } else if ((approximation == HERMITE) && (singularity == SKIPPING)) {
        for (int i = 0; i < tree_numnodes; i++)
            pc_comp_ms_modifiedF_hermite(tree, i, interpolationDegree, totalNumberInterpolationPoints,
                                         xS, yS, zS, qS, wS, xC, yC, zC, qC);

    } else if ((approximation == HERMITE) && (singularity == SUBTRACTION)) {
        for (int i = 0; i < tree_numnodes; i++)
            pc_comp_ms_modifiedF_hermite_SS(tree, i, interpolationDegree, totalNumberInterpolationPoints,
                                            xS, yS, zS, qS, wS, xC, yC, zC, qC, wC);

    } else {
        exit(1);
    }
    
    return;
}




void Clusters_Targets_Construct(struct Clusters **clusters_addr, const struct Particles *targets,
                                const struct Tree *tree, const struct RunParams *run_params)
{
    *clusters_addr = malloc(sizeof(struct Clusters));
    struct Clusters *clusters = *clusters_addr;

    SINGULARITY singularity = run_params->singularity;
    APPROXIMATION approximation = run_params->approximation;

    int tree_numnodes = tree->numnodes;
    int totalNumberTargetPoints = targets->num;

    int interpolationDegree = run_params->interp_degree;
    int interpDegreeLim = interpolationDegree + 1;
    int interpolationPointsPerCluster = run_params->interp_pts_per_cluster;

    int totalNumberInterpolationPoints  = tree_numnodes * interpolationPointsPerCluster;
    int totalNumberInterpolationCharges = tree_numnodes * run_params->interp_charges_per_cluster;
    int totalNumberInterpolationWeights = tree_numnodes * run_params->interp_weights_per_cluster;

    clusters->x = NULL;
    clusters->y = NULL;
    clusters->z = NULL;
    clusters->q = NULL;
    clusters->w = NULL;

    make_vector(clusters->x, totalNumberInterpolationPoints);
    make_vector(clusters->y, totalNumberInterpolationPoints);
    make_vector(clusters->z, totalNumberInterpolationPoints);
    make_vector(clusters->q, totalNumberInterpolationCharges);

    if (singularity == SUBTRACTION) {
        make_vector(clusters->w, totalNumberInterpolationWeights);
    }

    clusters->num         = totalNumberInterpolationPoints;
    clusters->num_charges = totalNumberInterpolationCharges;
    clusters->num_weights = totalNumberInterpolationWeights;

    double *xC = clusters->x;
    double *yC = clusters->y;
    double *zC = clusters->z;
    double *qC = clusters->q;
    double *wC = clusters->w;

#ifdef OPENACC_ENABLED
    #pragma acc enter data create(xC[0:totalNumberInterpolationPoints], yC[0:totalNumberInterpolationPoints], \
                                  zC[0:totalNumberInterpolationPoints], qC[0:totalNumberInterpolationCharges])

    #pragma acc kernels present(xC,yC,zC,qC)
    {
        for (int i=0;i<totalNumberInterpolationPoints;i++){
            xC[i]=0.0;
            yC[i]=0.0;
            zC[i]=0.0;
        }

        for (int i=0;i<totalNumberInterpolationCharges;i++){
            qC[i]=0.0;
        }
    }


    if (singularity == SUBTRACTION) {
        #pragma acc enter data create(wC[0:totalNumberInterpolationWeights])
        #pragma acc kernels present(wC)
        {
            for (int i=0;i<totalNumberInterpolationWeights;i++){
                wC[i]=0.0;
            }
        }
    }
#endif

    for (int i = 0; i < tree_numnodes; i++) {
        cp_comp_interp(tree, i, interpolationDegree, xC, yC, zC);
    }

#ifdef OPENACC_ENABLED
    #pragma acc wait
#endif

    return;
}




void Clusters_Alloc(struct Clusters **clusters_addr, int length, const struct RunParams *run_params)
{
    *clusters_addr = malloc(sizeof(struct Clusters));
    struct Clusters *clusters = *clusters_addr;

    APPROXIMATION approximation = run_params->approximation;
    SINGULARITY singularity = run_params->singularity;
    
    clusters->num = length;
    clusters->num_charges = length;
    clusters->num_weights = length;
    
    clusters->x = NULL;
    clusters->y = NULL;
    clusters->z = NULL;
    clusters->q = NULL;
    clusters->w = NULL;

    if (approximation == HERMITE)
        clusters->num_charges *= 8;

    if ((approximation == HERMITE) && (singularity == SUBTRACTION))
        clusters->num_weights *= 8;

    if (clusters->num > 0) {
        make_vector(clusters->x, clusters->num);
        make_vector(clusters->y, clusters->num);
        make_vector(clusters->z, clusters->num);
        make_vector(clusters->q, clusters->num_charges);
        if (singularity == SUBTRACTION) {
            make_vector(clusters->w, clusters->num_weights);
        }
    }

    return;
}   /* END of function allocate_cluster */




void Clusters_Free(struct Clusters **clusters_addr)
{
    struct Clusters *clusters = *clusters_addr;
    
    if (clusters != NULL) {
        if (clusters->x != NULL) free_vector(clusters->x);
        if (clusters->y != NULL) free_vector(clusters->y);
        if (clusters->z != NULL) free_vector(clusters->z);
        if (clusters->q != NULL) free_vector(clusters->q);
        if (clusters->w != NULL) free_vector(clusters->w);
        free(clusters);
    }
    
    clusters = NULL;

    return;
}   /* END of function Clusters_Free */




void Clusters_Free_Win(struct Clusters **clusters_addr)
{
    struct Clusters *clusters = *clusters_addr;

    if (clusters != NULL) {
        if (clusters->x != NULL) MPI_Free_mem(clusters->x);
        if (clusters->y != NULL) MPI_Free_mem(clusters->y);
        if (clusters->z != NULL) MPI_Free_mem(clusters->z);
        if (clusters->q != NULL) MPI_Free_mem(clusters->q);
        if (clusters->w != NULL) MPI_Free_mem(clusters->w);
        free(clusters);
    }

    clusters = NULL;

    return;
}   /* END of function Clusters_Free */




/************************************/
/***** LOCAL FUNCTIONS **************/
/************************************/

void pc_comp_ms_modifiedF_child_to_parent(const struct Tree *tree, int child_index, int parent_index, int interpolationDegree,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ)
{

    int interpDegreeLim = interpolationDegree + 1;
    int interpolationPointsPerCluster = interpDegreeLim * interpDegreeLim * interpDegreeLim;

    int child_startingIndexInClustersArray = child_index * interpolationPointsPerCluster;
    int parent_startingIndexInClustersArray = parent_index * interpolationPointsPerCluster;

    double *weights, *dj, *tt, *nodeX, *nodeY, *nodeZ, *modifiedF;
    int *exactIndX, *exactIndY, *exactIndZ;

    make_vector(weights,   interpDegreeLim);
    make_vector(dj,        interpDegreeLim);
    make_vector(tt,        interpDegreeLim);
    make_vector(nodeX,     interpDegreeLim);
    make_vector(nodeY,     interpDegreeLim);
    make_vector(nodeZ,     interpDegreeLim);
    make_vector(modifiedF, interpolationPointsPerCluster);
    make_vector(exactIndX, interpolationPointsPerCluster);
    make_vector(exactIndY, interpolationPointsPerCluster);
    make_vector(exactIndZ, interpolationPointsPerCluster);

    double x0 = tree->x_min[parent_index];
    double x1 = tree->x_max[parent_index];
    double y0 = tree->y_min[parent_index];
    double y1 = tree->y_max[parent_index];
    double z0 = tree->z_min[parent_index];
    double z1 = tree->z_max[parent_index];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(clusterX, clusterY, clusterZ, clusterQ) \
                       create(modifiedF[0:interpolationPointsPerCluster], exactIndX[0:interpolationPointsPerCluster], \
                              exactIndY[0:interpolationPointsPerCluster], exactIndZ[0:interpolationPointsPerCluster], \
                              nodeX[0:interpDegreeLim], nodeY[0:interpDegreeLim], \
                              nodeZ[0:interpDegreeLim], weights[0:interpDegreeLim], \
                              dj[0:interpDegreeLim], tt[0:interpDegreeLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpolationPointsPerCluster; j++) {
        modifiedF[j] = clusterQ[child_startingIndexInClustersArray + j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int i = 0; i < interpDegreeLim; i++) {
        tt[i] = cos(i * M_PI / interpolationDegree);
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpDegreeLim; j++) {
        dj[j] = 1.0;
        if (j == 0) dj[j] = 0.5;
        if (j == interpolationDegree) dj[j] = 0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpDegreeLim; j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

    // Compute modified f values
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interpolationPointsPerCluster; i++) { // loop through the source points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double sx = clusterX[child_startingIndexInClustersArray+i];
        double sy = clusterY[child_startingIndexInClustersArray+i];
        double sz = clusterZ[child_startingIndexInClustersArray+i];

#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) reduction(+:sumX) reduction(+:sumY) reduction(+:sumZ)
#endif
        for (int j = 0; j < (interpolationDegree+1); j++) {  // loop through the degree

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
        int k1 = j%(interpolationDegree+1);
        int kk = (j-k1)/(interpolationDegree+1);
        int k2 = kk%(interpolationDegree+1);
        kk = kk - k2;
        int k3 = kk / (interpolationDegree+1);

        double cz = nodeZ[k3];
        double w3 = weights[k3];

        double cy = nodeY[k2];
        double w2 = weights[k2];

        double cx = nodeX[k1];
        double w1 = weights[k1];

        // Fill cluster X, Y, and Z arrays
        clusterX[parent_startingIndexInClustersArray + j] = cx;
        clusterY[parent_startingIndexInClustersArray + j] = cy;
        clusterZ[parent_startingIndexInClustersArray + j] = cz;

        // Increment cluster Q array
        double temp = 0.0;
#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) reduction(+:temp)
#endif
        for (int i = 0; i < interpolationPointsPerCluster; i++) {  // loop over source points
            double sx = clusterX[child_startingIndexInClustersArray + i];
            double sy = clusterY[child_startingIndexInClustersArray + i];
            double sz = clusterZ[child_startingIndexInClustersArray + i];

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

        clusterQ[parent_startingIndexInClustersArray + j] += temp;

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


void pc_comp_ms_modifiedF_SS_child_to_parent(const struct Tree *tree, int child_index, int parent_index, int interpolationDegree,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW)
{

    int interpDegreeLim = interpolationDegree + 1;
    int interpolationPointsPerCluster = interpDegreeLim * interpDegreeLim * interpDegreeLim;

    int child_startingIndexInClustersArray = child_index * interpolationPointsPerCluster;
    int parent_startingIndexInClustersArray = parent_index * interpolationPointsPerCluster;

    double *weights, *dj, *tt, *nodeX, *nodeY, *nodeZ, *modifiedF, *modifiedF2;
    int *exactIndX, *exactIndY, *exactIndZ;

    make_vector(weights,   interpDegreeLim);
    make_vector(dj,        interpDegreeLim);
    make_vector(tt,        interpDegreeLim);
    make_vector(nodeX,     interpDegreeLim);
    make_vector(nodeY,     interpDegreeLim);
    make_vector(nodeZ,     interpDegreeLim);
    make_vector(modifiedF, interpolationPointsPerCluster);
    make_vector(modifiedF2, interpolationPointsPerCluster);
    make_vector(exactIndX, interpolationPointsPerCluster);
    make_vector(exactIndY, interpolationPointsPerCluster);
    make_vector(exactIndZ, interpolationPointsPerCluster);

    double x0 = tree->x_min[parent_index];
    double x1 = tree->x_max[parent_index];
    double y0 = tree->y_min[parent_index];
    double y1 = tree->y_max[parent_index];
    double z0 = tree->z_min[parent_index];
    double z1 = tree->z_max[parent_index];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(clusterX, clusterY, clusterZ, clusterQ, clusterW) \
                       create(modifiedF[0:interpolationPointsPerCluster], modifiedF2[0:interpolationPointsPerCluster], exactIndX[0:interpolationPointsPerCluster], \
                              exactIndY[0:interpolationPointsPerCluster], exactIndZ[0:interpolationPointsPerCluster], \
                              nodeX[0:interpDegreeLim], nodeY[0:interpDegreeLim], \
                              nodeZ[0:interpDegreeLim], weights[0:interpDegreeLim], \
                              dj[0:interpDegreeLim], tt[0:interpDegreeLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpolationPointsPerCluster; j++) {
        modifiedF[j] = clusterQ[child_startingIndexInClustersArray + j];
        modifiedF2[j] = clusterW[child_startingIndexInClustersArray + j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int i = 0; i < interpDegreeLim; i++) {
        tt[i] = cos(i * M_PI / interpolationDegree);
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpDegreeLim; j++) {
        dj[j] = 1.0;
        if (j == 0) dj[j] = 0.5;
        if (j == interpolationDegree) dj[j] = 0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpDegreeLim; j++) {
        weights[j] = ((j % 2 == 0)? 1 : -1) * dj[j];
    }

    // Compute modified f values
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < interpolationPointsPerCluster; i++) { // loop through the source points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double sx = clusterX[child_startingIndexInClustersArray+i];
        double sy = clusterY[child_startingIndexInClustersArray+i];
        double sz = clusterZ[child_startingIndexInClustersArray+i];

#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) reduction(+:sumX) reduction(+:sumY) reduction(+:sumZ)
#endif
        for (int j = 0; j < (interpolationDegree+1); j++) {  // loop through the degree

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
        modifiedF2[i] /= denominator;
    }


#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < interpolationPointsPerCluster; j++) {
        int k1 = j%(interpolationDegree+1);
        int kk = (j-k1)/(interpolationDegree+1);
        int k2 = kk%(interpolationDegree+1);
        kk = kk - k2;
        int k3 = kk / (interpolationDegree+1);

        double cz = nodeZ[k3];
        double w3 = weights[k3];

        double cy = nodeY[k2];
        double w2 = weights[k2];

        double cx = nodeX[k1];
        double w1 = weights[k1];

        // Fill cluster X, Y, and Z arrays
        clusterX[parent_startingIndexInClustersArray + j] = cx;
        clusterY[parent_startingIndexInClustersArray + j] = cy;
        clusterZ[parent_startingIndexInClustersArray + j] = cz;

        // Increment cluster Q array
        double temp = 0.0;
        double temp2 = 0.0;
#ifdef OPENACC_ENABLED
        #pragma acc loop vector(32) reduction(+:temp) reduction(+:temp2)
#endif
        for (int i = 0; i < interpolationPointsPerCluster; i++) {  // loop over source points
            double sx = clusterX[child_startingIndexInClustersArray + i];
            double sy = clusterY[child_startingIndexInClustersArray + i];
            double sz = clusterZ[child_startingIndexInClustersArray + i];

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
            temp2 += numerator * modifiedF2[i];

        }

        clusterQ[parent_startingIndexInClustersArray + j] += temp;
        clusterW[parent_startingIndexInClustersArray + j] += temp2;

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
    free_vector(modifiedF2);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}




void pc_comp_ms_modifiedF(const struct Tree *tree, int idx, int interpolationDegree,
        double *xS, double *yS, double *zS, double *qS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ)
{

    int interpDegreeLim = interpolationDegree + 1;
    int interpolationPointsPerCluster = interpDegreeLim * interpDegreeLim * interpDegreeLim;
    int sourcePointsInCluster = tree->iend[idx] - tree->ibeg[idx] + 1;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInSourcesArray = tree->ibeg[idx]-1;

    double *weights, *dj, *tt, *nodeX, *nodeY, *nodeZ, *modifiedF;
    int *exactIndX, *exactIndY, *exactIndZ;

    make_vector(weights,   interpDegreeLim);
    make_vector(dj,        interpDegreeLim);
    make_vector(tt,        interpDegreeLim);
    make_vector(nodeX,     interpDegreeLim);
    make_vector(nodeY,     interpDegreeLim);
    make_vector(nodeZ,     interpDegreeLim);
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
                              nodeX[0:interpDegreeLim], nodeY[0:interpDegreeLim], \
                              nodeZ[0:interpDegreeLim], weights[0:interpDegreeLim], \
                              dj[0:interpDegreeLim], tt[0:interpDegreeLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < sourcePointsInCluster; j++) {
        modifiedF[j] = qS[startingIndexInSourcesArray + j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int i = 0; i < interpDegreeLim; i++) {
        tt[i] = cos(i * M_PI / interpolationDegree);
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpDegreeLim; j++) {
        dj[j] = 1.0;
        if (j == 0) dj[j] = 0.5;
        if (j == interpolationDegree) dj[j] = 0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpDegreeLim; j++) {
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
        #pragma acc loop vector(32) reduction(+:sumX) reduction(+:sumY) reduction(+:sumZ)
#endif
        for (int j = 0; j < (interpolationDegree+1); j++) {  // loop through the degree

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
        int k1 = j%(interpolationDegree+1);
        int kk = (j-k1)/(interpolationDegree+1);
        int k2 = kk%(interpolationDegree+1);
        kk = kk - k2;
        int k3 = kk / (interpolationDegree+1);

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
        #pragma acc loop vector(32) reduction(+:temp)
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



void pc_comp_ms_modifiedF_SS(const struct Tree *tree, int idx, int interpolationDegree,
        double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW)
{
    int interpDegreeLim = interpolationDegree + 1;
    int pointsPerCluster =  interpDegreeLim * interpDegreeLim * interpDegreeLim;
    int pointsInNode = tree->iend[idx] - tree->ibeg[idx] + 1;
    int startingIndexInClusters = idx * pointsPerCluster;
    int startingIndexInSources = tree->ibeg[idx]-1;

    double *weights, *dj, *tt, *nodeX, *nodeY, *nodeZ, *modifiedF, *modifiedF2;
    int *exactIndX, *exactIndY, *exactIndZ;

    make_vector(weights,    interpDegreeLim);
    make_vector(dj,         interpDegreeLim);
    make_vector(tt,         interpDegreeLim);
    make_vector(nodeX,      interpDegreeLim);
    make_vector(nodeY,      interpDegreeLim);
    make_vector(nodeZ,      interpDegreeLim);
    make_vector(modifiedF,  pointsInNode);
    make_vector(modifiedF2, pointsInNode);
    make_vector(exactIndX,  pointsInNode);
    make_vector(exactIndY,  pointsInNode);
    make_vector(exactIndZ,  pointsInNode);

    double x0 = tree->x_min[idx];  // 1e-15 fails for large meshes, mysteriously.
    double x1 = tree->x_max[idx];
    double y0 = tree->y_min[idx];
    double y1 = tree->y_max[idx];
    double z0 = tree->z_min[idx];
    double z1 = tree->z_max[idx];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 3;
    #pragma acc kernels async(streamID) present(xS, yS, zS, qS, wS, \
                                                clusterX, clusterY, clusterZ, clusterQ, clusterW) \
                       create(modifiedF[0:pointsInNode], modifiedF2[0:pointsInNode], exactIndX[0:pointsInNode], \
                              exactIndY[0:pointsInNode], exactIndZ[0:pointsInNode], \
                              nodeX[0:interpDegreeLim], nodeY[0:interpDegreeLim], nodeZ[0:interpDegreeLim], \
                              weights[0:interpDegreeLim], dj[0:interpDegreeLim], tt[0:interpDegreeLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
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
    #pragma acc loop vector(32) independent
#endif
    for (int i = 0; i < interpDegreeLim; i++) {
        tt[i] = cos(i * M_PI / interpolationDegree);
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpDegreeLim; j++) {
        dj[j] = 1.0;
        if (j == 0) dj[j] = 0.5;
        if (j == interpolationDegree) dj[j] = 0.5;
    }

#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpDegreeLim; j++) {
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
        #pragma acc loop vector(32) reduction(+:sumX) reduction(+:sumY) reduction(+:sumZ)
#endif
        for (int j = 0; j < interpDegreeLim; j++) {  // loop through the degree

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
        int k1 = j % interpDegreeLim;
        int kk = (j-k1) / interpDegreeLim;
        int k2 = kk % interpDegreeLim;
        kk = kk - k2;
        int k3 = kk / interpDegreeLim;

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
        #pragma acc loop vector(32) reduction(+:temp) reduction(+:temp2)
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

    free_vector(weights);
    free_vector(dj);
    free_vector(tt);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);
    free_vector(modifiedF);
    free_vector(modifiedF2);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}




void pc_comp_ms_modifiedF_hermite(const struct Tree *tree, int idx, int interpolationDegree,
        int totalNumberInterpolationPoints, double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ)
{

    int interpDegreeLim = interpolationDegree + 1;
    int interpolationPointsPerCluster =  interpDegreeLim * interpDegreeLim * interpDegreeLim;
    int sourcePointsInCluster = tree->iend[idx] - tree->ibeg[idx] + 1;
    int startingIndexInSourcesArray = tree->ibeg[idx] - 1;

    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInClusterChargesArray = idx * interpolationPointsPerCluster * 8;

    double *dj, *tt, *ww, *wx, *wy, *wz;
    double *nodeX, *nodeY, *nodeZ, *modifiedF;
    int *exactIndX, *exactIndY, *exactIndZ;

    make_vector(dj,         interpDegreeLim);
    make_vector(tt,         interpDegreeLim);
    make_vector(ww,         interpDegreeLim);
    make_vector(wx,         interpDegreeLim);
    make_vector(wy,         interpDegreeLim);
    make_vector(wz,         interpDegreeLim);
    make_vector(nodeX,      interpDegreeLim);
    make_vector(nodeY,      interpDegreeLim);
    make_vector(nodeZ,      interpDegreeLim);
    make_vector(modifiedF,  sourcePointsInCluster);
    make_vector(exactIndX,  sourcePointsInCluster);
    make_vector(exactIndY,  sourcePointsInCluster);
    make_vector(exactIndZ,  sourcePointsInCluster);

    // Set the bounding box.
    double x0 = tree->x_min[idx];
    double x1 = tree->x_max[idx];
    double y0 = tree->y_min[idx];
    double y1 = tree->y_max[idx];
    double z0 = tree->z_min[idx];
    double z1 = tree->z_max[idx];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 3;
    #pragma acc kernels async(streamID) present(xS, yS, zS, qS, \
                                                clusterX, clusterY, clusterZ, clusterQ) \
                       create(modifiedF[0:sourcePointsInCluster], exactIndX[0:sourcePointsInCluster], \
                              exactIndY[0:sourcePointsInCluster], exactIndZ[0:sourcePointsInCluster], \
                              nodeX[0:interpDegreeLim], nodeY[0:interpDegreeLim], nodeZ[0:interpDegreeLim], \
                              dj[0:interpDegreeLim], tt[0:interpDegreeLim], ww[0:interpDegreeLim], \
                              wx[0:interpDegreeLim], wy[0:interpDegreeLim], wz[0:interpDegreeLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < sourcePointsInCluster; j++) {
        modifiedF[j] = qS[startingIndexInSourcesArray + j];
        exactIndX[j] = -1;
        exactIndY[j] = -1;
        exactIndZ[j] = -1;
    }

    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int i = 0; i < interpDegreeLim; i++) {
        double xx = i * M_PI / interpolationDegree;
        tt[i] =  cos(xx);
        ww[i] = -cos(xx) / (2 * sin(xx) * sin(xx));
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }
    ww[0] = 0.25 * (interpolationDegree*interpolationDegree/3.0 + 1.0/6.0);
    ww[interpolationDegree] = -ww[0];

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpDegreeLim; j++) {
        dj[j] = 1.0;
        wx[j] = -4.0 * ww[j] / (x1 - x0);
        wy[j] = -4.0 * ww[j] / (y1 - y0);
        wz[j] = -4.0 * ww[j] / (z1 - z0);
    }
    dj[0] = 0.25;
    dj[interpolationDegree] = 0.25;



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
        #pragma acc loop vector(32) reduction(+:sumX) reduction(+:sumY) reduction(+:sumZ)
#endif
        for (int j = 0; j < interpDegreeLim; j++) {  // loop through the degree

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
        int k1 = j % interpDegreeLim;
        int kk = (j-k1) / interpDegreeLim;
        int k2 = kk % interpDegreeLim;
        kk = kk - k2;
        int k3 = kk / interpDegreeLim;

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
        #pragma acc loop vector(32) reduction(+:temp0) reduction(+:temp1) reduction(+:temp2) \
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
                numerator0 *=  Ay;                     // aAa

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

        clusterQ[startingIndexInClusterChargesArray + 0 * interpolationPointsPerCluster + j] += temp0;
        clusterQ[startingIndexInClusterChargesArray + 1 * interpolationPointsPerCluster + j] += temp1;
        clusterQ[startingIndexInClusterChargesArray + 2 * interpolationPointsPerCluster + j] += temp2;
        clusterQ[startingIndexInClusterChargesArray + 3 * interpolationPointsPerCluster + j] += temp3;
        clusterQ[startingIndexInClusterChargesArray + 4 * interpolationPointsPerCluster + j] += temp4;
        clusterQ[startingIndexInClusterChargesArray + 5 * interpolationPointsPerCluster + j] += temp5;
        clusterQ[startingIndexInClusterChargesArray + 6 * interpolationPointsPerCluster + j] += temp6;
        clusterQ[startingIndexInClusterChargesArray + 7 * interpolationPointsPerCluster + j] += temp7;

    }

#ifdef OPENACC_ENABLED
    } // end acc kernels region
#endif

    free_vector(dj);
    free_vector(tt);
    free_vector(ww);
    free_vector(wx);
    free_vector(wy);
    free_vector(wz);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);
    free_vector(modifiedF);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}



void pc_comp_ms_modifiedF_hermite_SS(const struct Tree *tree, int idx, int interpolationDegree,
        int totalNumberInterpolationPoints, double *xS, double *yS, double *zS, double *qS, double *wS,
        double *clusterX, double *clusterY, double *clusterZ, double *clusterQ, double *clusterW)
{

    int interpDegreeLim = interpolationDegree + 1;
    int interpolationPointsPerCluster =  interpDegreeLim * interpDegreeLim * interpDegreeLim;
    int sourcePointsInCluster = tree->iend[idx] - tree->ibeg[idx] + 1;
    int startingIndexInSourcesArray = tree->ibeg[idx] - 1;

    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInClusterWeightsArray = idx * interpolationPointsPerCluster * 8;
    int startingIndexInClusterChargesArray = idx * interpolationPointsPerCluster * 8;

    double *dj, *tt, *ww, *wx, *wy, *wz;
    double *nodeX, *nodeY, *nodeZ, *modifiedF, *modifiedF2;
    int *exactIndX, *exactIndY, *exactIndZ;

    make_vector(dj,          interpDegreeLim);
    make_vector(tt,          interpDegreeLim);
    make_vector(ww,          interpDegreeLim);
    make_vector(wx,          interpDegreeLim);
    make_vector(wy,          interpDegreeLim);
    make_vector(wz,          interpDegreeLim);
    make_vector(nodeX,       interpDegreeLim);
    make_vector(nodeY,       interpDegreeLim);
    make_vector(nodeZ,       interpDegreeLim);
    make_vector(modifiedF,   sourcePointsInCluster);
    make_vector(modifiedF2,  sourcePointsInCluster);
    make_vector(exactIndX,   sourcePointsInCluster);
    make_vector(exactIndY,   sourcePointsInCluster);
    make_vector(exactIndZ,   sourcePointsInCluster);

    // Set the bounding box.
    double x0 = tree->x_min[idx];
    double x1 = tree->x_max[idx];
    double y0 = tree->y_min[idx];
    double y1 = tree->y_max[idx];
    double z0 = tree->z_min[idx];
    double z1 = tree->z_max[idx];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 3;
    #pragma acc kernels async(streamID) present(xS, yS, zS, qS, wS, \
                                                clusterX, clusterY, clusterZ, clusterQ) \
                       create(modifiedF[0:sourcePointsInCluster], modifiedF2[0:sourcePointsInCluster], \
                              exactIndX[0:sourcePointsInCluster], exactIndY[0:sourcePointsInCluster], \
                              exactIndZ[0:sourcePointsInCluster], \
                              nodeX[0:interpDegreeLim], nodeY[0:interpDegreeLim], nodeZ[0:interpDegreeLim], \
                              dj[0:interpDegreeLim], tt[0:interpDegreeLim], ww[0:interpDegreeLim], \
                              wx[0:interpDegreeLim], wy[0:interpDegreeLim], wz[0:interpDegreeLim])
    {
#endif

#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
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
    #pragma acc loop vector(32) independent
#endif
    for (int i = 0; i < interpDegreeLim; i++) {
        double xx = i * M_PI / interpolationDegree;
        tt[i] =  cos(xx);
        ww[i] = -cos(xx) / (2 * sin(xx) * sin(xx));
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }
    ww[0] = 0.25 * (interpolationDegree*interpolationDegree/3.0 + 1.0/6.0);
    ww[interpolationDegree] = -ww[0];

    // Compute weights
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpDegreeLim; j++) {
        dj[j] = 1.0;
        wx[j] = -4.0 * ww[j] / (x1 - x0);
        wy[j] = -4.0 * ww[j] / (y1 - y0);
        wz[j] = -4.0 * ww[j] / (z1 - z0);
    }
    dj[0] = 0.25;
    dj[interpolationDegree] = 0.25;


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
        #pragma acc loop vector(32) reduction(+:sumX) reduction(+:sumY) reduction(+:sumZ)
#endif
        for (int j = 0; j < interpDegreeLim; j++) {  // loop through the degree

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
        int k1 = j % interpDegreeLim;
        int kk = (j-k1) / interpDegreeLim;
        int k2 = kk % interpDegreeLim;
        kk = kk - k2;
        int k3 = kk / interpDegreeLim;

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
        #pragma acc loop vector(32) reduction(+:tempq0) reduction(+:tempq1) reduction(+:tempq2) \
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

         clusterQ[startingIndexInClusterChargesArray + 0 * interpolationPointsPerCluster + j] += tempq0;
         clusterQ[startingIndexInClusterChargesArray + 1 * interpolationPointsPerCluster + j] += tempq1;
         clusterQ[startingIndexInClusterChargesArray + 2 * interpolationPointsPerCluster + j] += tempq2;
         clusterQ[startingIndexInClusterChargesArray + 3 * interpolationPointsPerCluster + j] += tempq3;
         clusterQ[startingIndexInClusterChargesArray + 4 * interpolationPointsPerCluster + j] += tempq4;
         clusterQ[startingIndexInClusterChargesArray + 5 * interpolationPointsPerCluster + j] += tempq5;
         clusterQ[startingIndexInClusterChargesArray + 6 * interpolationPointsPerCluster + j] += tempq6;
         clusterQ[startingIndexInClusterChargesArray + 7 * interpolationPointsPerCluster + j] += tempq7;

         clusterW[startingIndexInClusterWeightsArray + 0 * interpolationPointsPerCluster + j] += tempw0;
         clusterW[startingIndexInClusterWeightsArray + 1 * interpolationPointsPerCluster + j] += tempw1;
         clusterW[startingIndexInClusterWeightsArray + 2 * interpolationPointsPerCluster + j] += tempw2;
         clusterW[startingIndexInClusterWeightsArray + 3 * interpolationPointsPerCluster + j] += tempw3;
         clusterW[startingIndexInClusterWeightsArray + 4 * interpolationPointsPerCluster + j] += tempw4;
         clusterW[startingIndexInClusterWeightsArray + 5 * interpolationPointsPerCluster + j] += tempw5;
         clusterW[startingIndexInClusterWeightsArray + 6 * interpolationPointsPerCluster + j] += tempw6;
         clusterW[startingIndexInClusterWeightsArray + 7 * interpolationPointsPerCluster + j] += tempw7;
    }

#ifdef OPENACC_ENABLED
    } // end acc kernels region
#endif

    free_vector(dj);
    free_vector(tt);
    free_vector(ww);
    free_vector(wx);
    free_vector(wy);
    free_vector(wz);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);
    free_vector(modifiedF);
    free_vector(modifiedF2);
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}



void cp_comp_interp(const struct Tree *tree, int idx, int interpolationDegree,
                    double *clusterX, double *clusterY, double *clusterZ)
{

    int interpDegreeLim = interpolationDegree + 1;
    int interpolationPointsPerCluster = interpDegreeLim * interpDegreeLim * interpDegreeLim;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;

    double *tt, *nodeX, *nodeY, *nodeZ;

    make_vector(tt,    interpDegreeLim);
    make_vector(nodeX, interpDegreeLim);
    make_vector(nodeY, interpDegreeLim);
    make_vector(nodeZ, interpDegreeLim);

    double x0 = tree->x_min[idx];
    double x1 = tree->x_max[idx];
    double y0 = tree->y_min[idx];
    double y1 = tree->y_max[idx];
    double z0 = tree->z_min[idx];
    double z1 = tree->z_max[idx];

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(clusterX, clusterY, clusterZ) \
                       create(nodeX[0:interpDegreeLim], nodeY[0:interpDegreeLim], \
                              nodeZ[0:interpDegreeLim], tt[0:interpDegreeLim])
    {
#endif


    //  Fill in arrays of unique x, y, and z coordinates for the interpolation points.
#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent 
#endif
    for (int i = 0; i < interpDegreeLim; i++) {
        tt[i] = cos(i * M_PI / interpolationDegree);
        nodeX[i] = x0 + (tt[i] + 1.0)/2.0 * (x1 - x0);
        nodeY[i] = y0 + (tt[i] + 1.0)/2.0 * (y1 - y0);
        nodeZ[i] = z0 + (tt[i] + 1.0)/2.0 * (z1 - z0);
    }


#ifdef OPENACC_ENABLED
    #pragma acc loop vector(32) independent
#endif
    for (int j = 0; j < interpolationPointsPerCluster; j++) {
        int k1 = j%(interpolationDegree+1);
        int kk = (j-k1)/(interpolationDegree+1);
        int k2 = kk%(interpolationDegree+1);
        kk = kk - k2;
        int k3 = kk / (interpolationDegree+1);

        // Fill cluster X, Y, and Z arrays
        clusterX[startingIndexInClustersArray + j] = nodeX[k1];
        clusterY[startingIndexInClustersArray + j] = nodeY[k2];
        clusterZ[startingIndexInClustersArray + j] = nodeZ[k3];
    }
#ifdef OPENACC_ENABLED
    } //end acc kernels region
#endif

    free_vector(tt);
    free_vector(nodeX);
    free_vector(nodeY);
    free_vector(nodeZ);

    return;
}
