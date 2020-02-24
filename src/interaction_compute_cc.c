
/*
 *Procedures for Cluster-Cluster Treecode
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


static void cc_comp_pot(struct tnode_array *tree_array, int idx, double *pointwisePotential, int interpolationOrder,
                        double *xT, double *yT, double *zT, double *qT,
                        double *clusterQ, double *clusterW);

//static void cp_comp_pot_SS(struct tnode_array *tree_array, int idx, int interpolationOrder,
//                      double *xT, double *yT, double *zT, double *qT,
//                      double *clusterQ, double *clusterW);

static void cc_comp_pot_hermite(struct tnode_array *tree_array, int idx, double *pointwisePotential, int interpolationOrder,
                        int totalNumberInterpolationPoints,
                        double *xT, double *yT, double *zT, double *qT,
                        double *clusterQ, double *clusterW);

//static void cp_comp_pot_hermite_SS(struct tnode_array *tree_array, int idx, int interpolationOrder,
//                      int totalNumberInterpolationPoints,
//                      double *xT, double *yT, double *zT, double *qT,
//                      double *clusterQ, double *clusterW);



void InteractionCompute_CC_1(struct tnode_array *source_tree_array, struct tnode_array *target_tree_array,
                             int *approx_inter_list, int *direct_inter_list,
                             double *source_x, double *source_y, double *source_z,
                             double *source_q, double *source_w,
                             double *target_x, double *target_y, double *target_z, double *target_q,
                             double *source_cluster_x, double *source_cluster_y, double *source_cluster_z,
                             double *source_cluster_q, double *source_cluster_w,
                             double *target_cluster_x, double *target_cluster_y, double *target_cluster_z,
                             double *target_cluster_q, double *target_cluster_w,
                             double *pointwisePotential, int interpolationOrder,
                             int numSources, int numTargets, int numSourceClusterPoints, int numTargetClusterPoints,
                             struct kernel *kernel, char *singularityHandling, char *approximationName)
{

    int source_tree_numnodes = source_tree_array->numnodes;
    int target_tree_numnodes = target_tree_array->numnodes;
    
    double *xS = source_x;
    double *yS = source_y;
    double *zS = source_z;
    double *qS = source_q;
    double *wS = source_w;

    double *xT = target_x;
    double *yT = target_y;
    double *zT = target_z;
    double *qT = target_q;

    double *xCS = source_cluster_x;
    double *yCS = source_cluster_y;
    double *zCS = source_cluster_z;
    double *qCS = source_cluster_q;
    double *wCS = source_cluster_w;
    
    double *xCT = target_cluster_x;
    double *yCT = target_cluster_y;
    double *zCT = target_cluster_z;
    double *qCT = target_cluster_q;
    double *wCT = target_cluster_w;
    
    int *source_tree_ibeg = source_tree_array->ibeg;
    int *source_tree_iend = source_tree_array->iend;
    int *source_tree_cluster_ind = source_tree_array->cluster_ind;
    
    int *target_tree_ibeg = target_tree_array->ibeg;
    int *target_tree_iend = target_tree_array->iend;
    int *target_tree_cluster_ind = target_tree_array->cluster_ind;
    
    int numSourceClusterCharges = numSourceClusterPoints;
    int numSourceClusterWeights = numSourceClusterPoints;
    
    int numTargetClusterCharges = numTargetClusterPoints;
    int numTargetClusterWeights = numTargetClusterPoints;

    if (strcmp(approximationName, "hermite") == 0) {
        numSourceClusterCharges = 8 * numSourceClusterPoints;
        numTargetClusterCharges = 8 * numTargetClusterPoints;
    }

    if ((strcmp(approximationName, "hermite") == 0) && (strcmp(singularityHandling, "subtraction") == 0)) {
        numSourceClusterWeights = 8 * numSourceClusterPoints;
        numTargetClusterWeights = 8 * numTargetClusterPoints;
    }
        
    //for (int i = 0; i < numTargets; i++) pointwisePotential[i] = 0.0;

#ifdef OPENACC_ENABLED
    #pragma acc data copyin(xS[0:numSources], yS[0:numSources], zS[0:numSources], \
                            qS[0:numSources], wS[0:numSources], \
                            xT[0:numTargets], yT[0:numTargets], zT[0:numTargets], \
                            qT[0:numTargets], \
                            xCS[0:numSourceClusterPoints], \
                            yCS[0:numSourceClusterPoints], \
                            zCS[0:numSourceClusterPoints], \
                            qCS[0:numSourceClusterCharges], \
                            wCS[0:numSourceClusterWeights], \
                            xCT[0:numTargetClusterPoints], \
                            yCT[0:numTargetClusterPoints], \
                            zCT[0:numTargetClusterPoints], \
                       copy(qCT[0:numTargetClusterCharges], \
                            wCT[0:numTargetClusterWeights], \
                            pointwisePotential[0:numTargets])
#endif
    {

    int numInterpPoints = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);

    for (int i = 0; i < target_tree_numnodes; i++) {
        int target_ibeg = target_tree_ibeg[i];
        int target_iend = target_tree_iend[i];
        int target_cluster_ind = target_tree_cluster_ind[i];
        int target_cluster_start = numInterpPoints * target_cluster_ind;

        int numTargets = target_iend - target_ibeg + 1;
        int targetStart =  target_ibeg - 1;
        
        int numClusterApprox = target_tree_numApprox[i];
        int numClusterDirect = target_tree_numDirect[i];


/**********************************************************/
/************** POTENTIAL FROM APPROX *********************/
/**********************************************************/

        for (int j = 0; j < numClusterApprox; j++) {
            int source_node_index = approx_inter_list[i][j];
            int source_cluster_start = numInterpPoints * source_cluster_ind[source_node_index];
            int streamID = j%3;


    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/
            if (strcmp(kernel->name, "coulomb") == 0) {

                if (strcmp(approximationName, "lagrange") == 0) {

                    if (strcmp(singularityHandling, "skipping") == 0) {
            
                        CP_coulombApproximationLagrange(numInterpPoints, numInterpPoints,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
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

                        CP_coulombApproximationHermite(numInterpPoints, numInterpPoints,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
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

            int source_node_index = direct_inter_list[i][j];
            int source_ibeg = source_tree_ibeg[source_node_index];
            int source_iend = source_tree_iend[source_node_index];
            
            int numSources = source_iend - source_ibeg + 1;
            int sourceStart =  source_ibeg - 1;
            int streamID = j%3;

    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/

            if (strcmp(kernel->name, "coulomb") == 0) {

                if (strcmp(singularityHandling, "skipping") == 0) {

                    coulombDirect(numTargets, numSources,
                            targetStart, sourceStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            kernel, pointwisePotential, streamID);

                } else if (strcmp(singularityHandling, "subtraction") == 0) {

                    coulombSingularitySubtractionDirect(numTargets, numSources,
                            targetStart, sourceStart,
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

                    yukawaDirect(numTargets, numSources,
                            targetStart, sourceStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            kernel, pointwisePotential, streamID);

                } else if (strcmp(singularityHandling, "subtraction") == 0) {

                    yukawaSingularitySubtractionDirect(numTargets, numSources,
                            targetStart, sourceStart,
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
