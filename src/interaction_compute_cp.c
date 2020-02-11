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
#include "struct_nodes.h"
#include "struct_particles.h"
#include "struct_kernel.h"

#include "kernels/coulomb.h"
#include "kernels/cp_coulomb.h"
#include "kernels/yukawa.h"
#include "kernels/coulomb_singularity_subtraction.h"
#include "kernels/yukawa_singularity_subtraction.h"
#include "kernels/regularized-coulomb.h"
#include "kernels/regularized-yukawa.h"

#include "interaction_compute.h"


static void cp_comp_pot(struct tnode_array *tree_array, int idx, double *pointwisePotential, int interpolationOrder,
                        double *xT, double *yT, double *zT, double *qT,
                        double *clusterQ, double *clusterW);

//static void cp_comp_pot_SS(struct tnode_array *tree_array, int idx, int interpolationOrder,
//                      double *xT, double *yT, double *zT, double *qT,
//                      double *clusterQ, double *clusterW);

static void cp_comp_pot_hermite(struct tnode_array *tree_array, int idx, double *pointwisePotential, int interpolationOrder,
                        int totalNumberInterpolationPoints,
                        double *xT, double *yT, double *zT, double *qT,
                        double *clusterQ, double *clusterW);

//static void cp_comp_pot_hermite_SS(struct tnode_array *tree_array, int idx, int interpolationOrder,
//                      int totalNumberInterpolationPoints,
//                      double *xT, double *yT, double *zT, double *qT,
//                      double *clusterQ, double *clusterW);



void Interaction_Compute_CP1(struct tnode_array *tree_array, struct tnode_array *batches,
                             int *tree_inter_list, int *direct_inter_list,
                             double *source_x, double *source_y, double *source_z,
                             double *source_q, double *source_w,
                             double *target_x, double *target_y, double *target_z, double *target_q,
                             double *cluster_x, double *cluster_y, double *cluster_z,
                             double *cluster_q, double *cluster_w,
                             double *pointwisePotential, int interpolationOrder,
                             int numSources, int numTargets, int totalNumberOfInterpolationPoints,
                             int batch_approx_offset, int batch_direct_offset,
                             struct kernel *kernel, char *singularityHandling,
                             char *approximationName)
{

    int tree_numnodes = tree_array->numnodes;
    int batch_numnodes = batches->numnodes;
    
    double *xS = source_x;
    double *yS = source_y;
    double *zS = source_z;
    double *qS = source_q;
    double *wS = source_w;

    double *xT = target_x;
    double *yT = target_y;
    double *zT = target_z;
    double *qT = target_q;

    double *xC = cluster_x;
    double *yC = cluster_y;
    double *zC = cluster_z;
    double *qC = cluster_q;
    double *wC = cluster_w;

    int *ibegs = tree_array->ibeg;
    int *iends = tree_array->iend;
    int *clusterInd = tree_array->cluster_ind;

    int numberOfClusterCharges = totalNumberOfInterpolationPoints;
    int numberOfClusterWeights = totalNumberOfInterpolationPoints;

    if (strcmp(approximationName, "hermite") == 0)
        numberOfClusterCharges = 8 * totalNumberOfInterpolationPoints;

    if ((strcmp(approximationName, "hermite") == 0) && (strcmp(singularityHandling, "subtraction") == 0))
        numberOfClusterWeights = 8 * totalNumberOfInterpolationPoints;
        
    //for (int i = 0; i < numTargets; i++) pointwisePotential[i] = 0.0;

#ifdef OPENACC_ENABLED
    #pragma acc data copyin(xS[0:numSources], yS[0:numSources], zS[0:numSources], \
                            qS[0:numSources], wS[0:numSources], \
                            xT[0:numTargets], yT[0:numTargets], zT[0:numTargets], \
                            qT[0:numTargets], \
                            xC[0:totalNumberOfInterpolationPoints], \
                            yC[0:totalNumberOfInterpolationPoints], \
                            zC[0:totalNumberOfInterpolationPoints], \
                            tree_inter_list[0:batch_approx_offset*batch_numnodes], \
                            direct_inter_list[0:batch_direct_offset*batch_numnodes], \
                            ibegs[0:tree_numnodes], iends[0:tree_numnodes]) \
                       copy(qC[0:numberOfClusterCharges], \
                            wC[0:numberOfClusterWeights], \
                            pointwisePotential[0:numTargets])
#endif
    {

    int numberOfInterpolationPoints = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);

    for (int i = 0; i < batches->numnodes; i++) {
        int batch_ibeg = batches->ibeg[i];
        int batch_iend = batches->iend[i];
        int numberOfClusterApproximations = batches->numApprox[i];
        int numberOfDirectSums = batches->numDirect[i];

        int numberOfSources = batch_iend - batch_ibeg + 1;
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
            if (strcmp(kernel->name, "coulomb") == 0) {

                if (strcmp(approximationName, "lagrange") == 0) {

                    if (strcmp(singularityHandling, "skipping") == 0) {
            
                        CP_coulombApproximationLagrange(numberOfSources,
                            numberOfInterpolationPoints, batchStart, clusterStart,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            kernel, streamID);

                    } else if (strcmp(singularityHandling, "subtraction") == 0) {

                        printf("Not yet implemented!\n");
                        exit(1);

                    } else {
                        printf("Invalid choice of singularityHandling. Exiting. \n");
                        exit(1);
                    }

                } else if (strcmp(approximationName, "hermite") == 0) {

                    if (strcmp(singularityHandling, "skipping") == 0) {

                        CP_coulombApproximationHermite(numberOfSources,
                            numberOfInterpolationPoints, batchStart, clusterStart,
                            totalNumberOfInterpolationPoints,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            kernel, streamID);

                    } else if (strcmp(singularityHandling, "subtraction") == 0) {

                        printf("Not yet implemented!\n");
                        exit(1);

                    } else {
                        printf("Invalid choice of singularityHandling. Exiting. \n");
                        exit(1);
                    }


                }else{
                    printf("Invalid approximationName.  Was set to %s\n", approximationName);
                    exit(1);
                }

    /***********************************************/
    /***************** Yukawa **********************/
    /***********************************************/

            } else if (strcmp(kernel->name, "yukawa") == 0) {

                if (strcmp(approximationName, "lagrange") == 0) {

                    if (strcmp(singularityHandling, "skipping") == 0) {

                        printf("Not yet implemented!\n");
                        exit(1);

                    } else if (strcmp(singularityHandling, "subtraction") == 0) {
                        
                        printf("Not yet implemented!\n");
                        exit(1);

                    } else {
                        printf("Invalid choice of singularityHandling. Exiting. \n");
                        exit(1);
                    }

                } else if (strcmp(approximationName, "hermite") == 0) {

                    if (strcmp(singularityHandling, "skipping") == 0) {

                        printf("Not yet implemented!\n");
                        exit(1);

                    } else if (strcmp(singularityHandling, "subtraction") == 0) {

                        printf("Not yet implemented!\n");
                        exit(1);

                    } else {
                        printf("Invalid choice of singularityHandling. Exiting. \n");
                        exit(1);
                    }

                } else {
                    printf("Invalid approximationName.\n");
                    exit(1);
                }

            } else {
                printf("Invalid kernel->name. Exiting.\n");
                exit(1);
            }

        } // end loop over cluster approximations



/**********************************************************/
/************** POTENTIAL FROM DIRECT *********************/
/**********************************************************/

        for (int j = 0; j < numberOfDirectSums; j++) {

            int node_index = direct_inter_list[i * batch_direct_offset + j];
            int target_start = ibegs[node_index] - 1;
            int target_end = iends[node_index];
            int number_of_targets_in_cluster = target_end-target_start;
            int streamID = j%3;

    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/

            if (strcmp(kernel->name, "coulomb") == 0) {

                if (strcmp(singularityHandling, "skipping") == 0) {

                    coulombDirect(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            kernel, pointwisePotential, streamID);

                } else if (strcmp(singularityHandling, "subtraction") == 0) {

                    coulombSingularitySubtractionDirect(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            kernel, pointwisePotential, streamID);

                }else {
                    printf("Invalid choice of singularityHandling. Exiting. \n");
                    exit(1);
                }

    /***********************************************/
    /***************** Yukawa **********************/
    /***********************************************/

            } else if (strcmp(kernel->name, "yukawa") == 0) {

                if (strcmp(singularityHandling, "skipping") == 0) {

                    yukawaDirect(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            kernel, pointwisePotential, streamID);

                } else if (strcmp(singularityHandling, "subtraction") == 0) {

                    yukawaSingularitySubtractionDirect(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            kernel, pointwisePotential, streamID);

                } else {
                    printf("Invalid choice of singularityHandling. Exiting. \n");
                    exit(1);
                }

            } else {
                printf("Invalid kernel->name. Exiting.\n");
                exit(1);
            }

        } // end loop over number of direct interactions

    } // end loop over target batches
#ifdef OPENACC_ENABLED
        #pragma acc wait
#endif
    } // end acc data region

    return;

} /* END of function pc_treecode */


void Interaction_Compute_CP2(struct tnode_array *tree_array,
                             double *target_x, double *target_y, double *target_z, double *target_q,
                             double *cluster_x, double *cluster_y, double *cluster_z,
                             double *cluster_q, double *cluster_w,
                             double *pointwisePotential, int interpolationOrder,
                             int numTargets, int totalNumberInterpolationPoints,
                             int totalNumberInterpolationCharges, int totalNumberInterpolationWeights,
                             char *singularityHandling, char *approximationName)
{

    int interpOrderLim = interpolationOrder+1;
    int tree_numnodes = tree_array->numnodes;

#ifdef OPENACC_ENABLED
    #pragma acc data copyin(tt[0:interpOrderLim], ww[0:interpOrderLim], \
                            target_x[0:numTargets], target_y[0:numTargets], \
                            target_z[0:numTargets], target_q[0:numTargets], \
                            cluster_q[0:totalNumberInterpolationCharges], \
                            cluster_w[0:totalNumberInterpolationWeights]) \
                       copy(pointwisePotential[0:numTargets])
    {
#endif

    if ((strcmp(approximationName, "lagrange") == 0) && (strcmp(singularityHandling, "skipping") == 0)) {
        for (int i = 0; i < tree_numnodes; i++)
            cp_comp_pot(tree_array, i, pointwisePotential, interpolationOrder,
                        target_x, target_y, target_z, target_q, cluster_q, cluster_w);

    } else if ((strcmp(approximationName, "lagrange") == 0) && (strcmp(singularityHandling, "subtraction") == 0)) {
//        for (int i = 0; i < tree_numnodes; i++)
//            cp_comp_pot_SS(tree_array, i, pointwisePotential interpolationOrder,
//                       target_x, target_y, target_z, target_q, cluster_q, cluster_w);

    } else if ((strcmp(approximationName, "hermite") == 0) && (strcmp(singularityHandling, "skipping") == 0)) {
        for (int i = 0; i < tree_numnodes; i++)
            cp_comp_pot_hermite(tree_array, i, pointwisePotential, interpolationOrder, totalNumberInterpolationPoints,
                        target_x, target_y, target_z, target_q, cluster_q, cluster_w);

    } else if ((strcmp(approximationName, "hermite") == 0) && (strcmp(singularityHandling, "subtraction") == 0)) {
//        for (int i = 0; i < tree_numnodes; i++)
//            cp_comp_pot_hermite_SS(tree_array, i, pointwisePotential, interpolationOrder, totalNumberInterpolationPoints,
//                       target_x, target_y, target_z, target_q, cluster_q, cluster_w);

    } else {
        exit(1);
    }
        
#ifdef OPENACC_ENABLED
    #pragma acc wait
    } // end ACC DATA REGION
#endif
    
    return;
}


void cp_comp_pot(struct tnode_array *tree_array, int idx, double *pointwisePotential, int interpolationOrder,
        double *target_x, double *target_y, double *target_z, double *target_q,
        double *cluster_q, double *cluster_w)
{
    int interpOrderLim = interpolationOrder + 1;
    int interpolationPointsPerCluster = interpOrderLim * interpOrderLim * interpOrderLim;
    int targetPointsInCluster = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInTargetsArray = tree_array->ibeg[idx]-1;
    
    double nodeX[interpOrderLim], nodeY[interpOrderLim], nodeZ[interpOrderLim];
    double weights[interpOrderLim], dj[interpOrderLim];
    int *exactIndX, *exactIndY, *exactIndZ;
    
    double x0 = tree_array->x_min[idx];
    double x1 = tree_array->x_max[idx];
    double y0 = tree_array->y_min[idx];
    double y1 = tree_array->y_max[idx];
    double z0 = tree_array->z_min[idx];
    double z1 = tree_array->z_max[idx];
    
    make_vector(exactIndX, targetPointsInCluster);
    make_vector(exactIndY, targetPointsInCluster);
    make_vector(exactIndZ, targetPointsInCluster);

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(target_x, target_y, target_z, target_q, cluster_q, tt) \
        create(exactIndX[0:targetPointsInCluster], exactIndY[0:targetPointsInCluster], exactIndZ[0:targetPointsInCluster], \
               nodeX[0:interpOrderLim], nodeY[0:interpOrderLim], nodeZ[0:interpOrderLim], \
               weights[0:interpOrderLim], dj[0:interpOrderLim])
    {
#endif
    
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < targetPointsInCluster; j++) {
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
    for (int j = 0; j < interpOrderLim; j++){
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

#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < targetPointsInCluster; i++) { // loop through the target points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double tx = target_x[startingIndexInTargetsArray+i];
        double ty = target_y[startingIndexInTargetsArray+i];
        double tz = target_z[startingIndexInTargetsArray+i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (int j = 0; j < interpOrderLim; j++) {  // loop through the degree

            double cx = tx - nodeX[j];
            double cy = ty - nodeY[j];
            double cz = tz - nodeZ[j];

            if (fabs(cx)<DBL_MIN) exactIndX[i] = j;
            if (fabs(cy)<DBL_MIN) exactIndY[i] = j;
            if (fabs(cz)<DBL_MIN) exactIndZ[i] = j;

            // Increment the sums
            double w = weights[j];
            sumX += w / cx;
            sumY += w / cy;
            sumZ += w / cz;

        }

        double denominator = 1.0;
        if (exactIndX[i]==-1) denominator /= sumX;
        if (exactIndY[i]==-1) denominator /= sumY;
        if (exactIndZ[i]==-1) denominator /= sumZ;
        
        double temp = 0.0;
        
#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp)
#endif
        for (int j = 0; j < interpolationPointsPerCluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point

            int k1 = j%interpOrderLim;
            int kk = (j-k1)/interpOrderLim;
            int k2 = kk%interpOrderLim;
            kk = kk - k2;
            int k3 = kk / interpOrderLim;

            double w3 = weights[k3];
            double w2 = weights[k2];
            double w1 = weights[k1];
            
            double cx = nodeX[k1];
            double cy = nodeY[k2];
            double cz = nodeZ[k3];
            double cq = cluster_q[startingIndexInClustersArray + j];
        
            double numerator = 1.0;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
            if (exactIndX[i] == -1) {
                numerator *= w1 / (tx - cx);
            } else {
                if (exactIndX[i] != k1) numerator *= 0;
            }

            if (exactIndY[i] == -1) {
                numerator *= w2 / (ty - cy);
            } else {
                if (exactIndY[i] != k2) numerator *= 0;
            }

            if (exactIndZ[i] == -1) {
                numerator *= w3 / (tz - cz);
            } else {
                if (exactIndZ[i] != k3) numerator *= 0;
            }

            temp += numerator * denominator * cq;
        }
        
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        pointwisePotential[i + startingIndexInTargetsArray] += temp;
    }
#ifdef OPENACC_ENABLED
    } //end ACC kernels
#endif
    
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}


void cp_comp_pot_hermite(struct tnode_array *tree_array, int idx, double *pointwisePotential, int interpolationOrder,
        int totalNumberInterpolationPoints, double *target_x, double *target_y, double *target_z, double *target_q,
        double *cluster_q, double *cluster_w)
{
    int interpOrderLim = interpolationOrder + 1;
    int interpolationPointsPerCluster = interpOrderLim * interpOrderLim * interpOrderLim;
    int targetPointsInCluster = tree_array->iend[idx] - tree_array->ibeg[idx] + 1;
    int startingIndexInClustersArray = idx * interpolationPointsPerCluster;
    int startingIndexInTargetsArray = tree_array->ibeg[idx]-1;
    
    double nodeX[interpOrderLim], nodeY[interpOrderLim], nodeZ[interpOrderLim];
    double wx[interpOrderLim], wy[interpOrderLim], wz[interpOrderLim], dj[interpOrderLim];
    int *exactIndX, *exactIndY, *exactIndZ;
    
    double *cluster_q_     = &cluster_q[8*startingIndexInClustersArray + 0*interpolationPointsPerCluster];
    double *cluster_q_dx   = &cluster_q[8*startingIndexInClustersArray + 1*interpolationPointsPerCluster];
    double *cluster_q_dy   = &cluster_q[8*startingIndexInClustersArray + 2*interpolationPointsPerCluster];
    double *cluster_q_dz   = &cluster_q[8*startingIndexInClustersArray + 3*interpolationPointsPerCluster];
    double *cluster_q_dxy  = &cluster_q[8*startingIndexInClustersArray + 4*interpolationPointsPerCluster];
    double *cluster_q_dyz  = &cluster_q[8*startingIndexInClustersArray + 5*interpolationPointsPerCluster];
    double *cluster_q_dxz  = &cluster_q[8*startingIndexInClustersArray + 6*interpolationPointsPerCluster];
    double *cluster_q_dxyz = &cluster_q[8*startingIndexInClustersArray + 7*interpolationPointsPerCluster];
    
    double x0 = tree_array->x_min[idx];
    double x1 = tree_array->x_max[idx];
    double y0 = tree_array->y_min[idx];
    double y1 = tree_array->y_max[idx];
    double z0 = tree_array->z_min[idx];
    double z1 = tree_array->z_max[idx];
    
    make_vector(exactIndX, targetPointsInCluster);
    make_vector(exactIndY, targetPointsInCluster);
    make_vector(exactIndZ, targetPointsInCluster);

#ifdef OPENACC_ENABLED
    int streamID = rand() % 4;
    #pragma acc kernels async(streamID) present(target_x, target_y, target_z, target_q, \
                                        cluster_q_, cluster_q_dx, cluster_q_dy, cluster_q_dz, \
                                        cluster_q_dxy, cluster_q_dyz, cluster_q_dxz, \
                                        cluster_q_dxyz, tt, ww) \
        create(exactIndX[0:targetPointsInCluster], exactIndY[0:targetPointsInCluster], exactIndZ[0:targetPointsInCluster], \
               nodeX[0:interpOrderLim], nodeY[0:interpOrderLim], nodeZ[0:interpOrderLim], \
               wx[0:interpOrderLim], wy[0:interpOrderLim], wz[0:interpOrderLim], dj[0:interpOrderLim])
    {
#endif
    
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int j = 0; j < targetPointsInCluster; j++) {
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
    for (int j = 0; j < interpOrderLim; j++){
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
    for (int i = 0; i < targetPointsInCluster; i++) { // loop through the target points

        double sumX = 0.0;
        double sumY = 0.0;
        double sumZ = 0.0;

        double tx = target_x[startingIndexInTargetsArray+i];
        double ty = target_y[startingIndexInTargetsArray+i];
        double tz = target_z[startingIndexInTargetsArray+i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent
#endif
        for (int j = 0; j < interpOrderLim; j++) {  // loop through the degree

            double cx =  tx - nodeX[j];
            double cy =  ty - nodeY[j];
            double cz =  tz - nodeZ[j];

            if (fabs(cx)<DBL_MIN) exactIndX[i] = j;
            if (fabs(cy)<DBL_MIN) exactIndY[i] = j;
            if (fabs(cz)<DBL_MIN) exactIndZ[i] = j;

            // Increment the sums
            sumX += dj[j] / (cx*cx) + wx[j] / cx;
            sumY += dj[j] / (cy*cy) + wy[j] / cy;
            sumZ += dj[j] / (cz*cz) + wz[j] / cz;

        }

        double denominator = 1.0;
        if (exactIndX[i]==-1) denominator /= sumX;
        if (exactIndY[i]==-1) denominator /= sumY;
        if (exactIndZ[i]==-1) denominator /= sumZ;
        
        double temp = 0.0;
        
#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp)
#endif
        for (int j = 0; j < interpolationPointsPerCluster; j++) { // loop over interpolation points, set (cx,cy,cz) for this point

            int k1 = j%interpOrderLim;
            int kk = (j-k1)/interpOrderLim;
            int k2 = kk%interpOrderLim;
            kk = kk - k2;
            int k3 = kk / interpOrderLim;
            
            double dx =  tx - nodeX[k1];
            double dy =  tx - nodeY[k2];
            double dz =  tx - nodeZ[k3];
            
            double cq     = cluster_q_[j];
            double cqdx   = cluster_q_dx[j];
            double cqdy   = cluster_q_dy[j];
            double cqdz   = cluster_q_dz[j];
            double cqdxy  = cluster_q_dxy[j];
            double cqdyz  = cluster_q_dyz[j];
            double cqdxz  = cluster_q_dxz[j];
            double cqdxyz = cluster_q_dxyz[j];
        
            double numerator0 = 1.0, numerator1 = 1.0, numerator2 = 1.0, numerator3 = 1.0;
            double numerator4 = 1.0, numerator5 = 1.0, numerator6 = 1.0, numerator7 = 1.0;

            double Ax = dj[k1] / (dx*dx) + wx[k1] / dx;
            double Ay = dj[k2] / (dy*dy) + wy[k2] / dy;
            double Az = dj[k3] / (dz*dz) + wz[k3] / dz;
            double Bx = dj[k1] / dx;
            double By = dj[k2] / dy;
            double Bz = dj[k3] / dz;

            // If exactInd[i] == -1, then no issues.
            // If exactInd[i] != -1, then we want to zero out terms EXCEPT when exactInd=k1.
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

            temp += denominator * (numerator0 * cq     +  numerator1 * cqdx   +  numerator2 * cqdy
                                +  numerator3 * cqdz   +  numerator4 * cqdxy  +  numerator5 * cqdyz
                                +  numerator6 * cqdxz  +  numerator7 * cqdxyz);
        }
        
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        pointwisePotential[i + startingIndexInTargetsArray] += temp;
    }
#ifdef OPENACC_ENABLED
    } //end ACC kernels
#endif
    
    free_vector(exactIndX);
    free_vector(exactIndY);
    free_vector(exactIndZ);

    return;
}
