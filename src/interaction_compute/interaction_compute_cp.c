/*
 *Procedures for Cluster-Particle Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "../utilities/array.h"

#include "../tree/struct_nodes.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"

#include "../kernels/coulomb/coulomb.h"
#include "../kernels/yukawa/yukawa.h"
#include "../kernels/regularized-coulomb/regularized-coulomb.h"
#include "../kernels/regularized-yukawa/regularized-yukawa.h"

#include "interaction_compute.h"


void InteractionCompute_CP(struct tnode_array *tree_array, struct tnode_array *batches,
                           int **approx_inter_list, int **direct_inter_list,
                           double *source_x, double *source_y, double *source_z,
                           double *source_q, double *source_w,
                           double *target_x, double *target_y, double *target_z, double *target_q,
                           double *cluster_x, double *cluster_y, double *cluster_z,
                           double *cluster_q, double *cluster_w,
                           double *pointwisePotential,
                           int numSources, int numTargets, int totalNumberOfInterpolationPoints,
                           struct RunParams *run_params)
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

    if (run_params->approximation == HERMITE)
        numberOfClusterCharges = 8 * totalNumberOfInterpolationPoints;

    if ((run_params->approximation == HERMITE) && (run_params->singularity == SUBTRACTION))
        numberOfClusterWeights = 8 * totalNumberOfInterpolationPoints;
        

#ifdef OPENACC_ENABLED
    #pragma acc data copyin(xS[0:numSources], yS[0:numSources], zS[0:numSources], \
                            qS[0:numSources], wS[0:numSources], \
                            xT[0:numTargets], yT[0:numTargets], zT[0:numTargets], \
                            qT[0:numTargets], \
                            xC[0:totalNumberOfInterpolationPoints], \
                            yC[0:totalNumberOfInterpolationPoints], \
                            zC[0:totalNumberOfInterpolationPoints]) \
                       copy(qC[0:numberOfClusterCharges], \
                            wC[0:numberOfClusterWeights], \
                            pointwisePotential[0:numTargets])
#endif
    {

    int numberOfInterpolationPoints = run_params->interp_pts_per_cluster;

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

            int node_index = approx_inter_list[i][j];

            int clusterStart = numberOfInterpolationPoints*clusterInd[node_index];
            int streamID = j%3;


    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/
            if (run_params->kernel == COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {
            
                        K_Coulomb_CP_Lagrange(numberOfSources,
                            numberOfInterpolationPoints, batchStart, clusterStart,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP COULOMB SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Coulomb_CP_Hermite(numberOfSources,
                            numberOfInterpolationPoints, batchStart, clusterStart,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP COULOMB SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else {
                    printf("**ERROR** INVALID CHOICE OF APPROXIMATION. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /***************** Yukawa **********************/
    /***********************************************/

            } else if (run_params->kernel == YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Yukawa_CP_Lagrange(numberOfSources,
                            numberOfInterpolationPoints, batchStart, clusterStart,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP YUKAWA SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Yukawa_CP_Hermite(numberOfSources,
                            numberOfInterpolationPoints, batchStart, clusterStart,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP YUKAWA SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else {
                    printf("**ERROR** INVALID CHOICE OF APPROXIMATION. EXITING. \n");
                    exit(1);
                }

    /***************************************/
    /********* Regularized Coulomb *********/
    /***************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_CP_Lagrange(numberOfSources,
                            numberOfInterpolationPoints, batchStart, clusterStart,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP REGULARIZED COULOMB SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_CP_Hermite(numberOfSources,
                            numberOfInterpolationPoints, batchStart, clusterStart,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP REGULARIZED COULOMB SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else {
                    printf("**ERROR** INVALID CHOICE OF APPROXIMATION. EXITING. \n");
                    exit(1);
                }

    /***************************************/
    /********* Regularized Yukawa **********/
    /***************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedYukawa_CP_Lagrange(numberOfSources,
                            numberOfInterpolationPoints, batchStart, clusterStart,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP REGULARIZED YUKAWA SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        printf("**ERROR** NOT SET UP FOR CP REGULARIZED YUKAWA HERMITE.  EXITING.\n");
                        exit(1);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP REGULARIZED COULOMB SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else {
                    printf("**ERROR** INVALID CHOICE OF APPROXIMATION. EXITING. \n");
                    exit(1);
                }

            } else {
                printf("**ERROR** INVALID KERNEL. EXITING.\n");
                exit(1);
            }

        } // end loop over cluster approximations



/**********************************************************/
/************** POTENTIAL FROM DIRECT *********************/
/**********************************************************/

        for (int j = 0; j < numberOfDirectSums; j++) {

            int node_index = direct_inter_list[i][j];

            int target_start = ibegs[node_index] - 1;
            int target_end = iends[node_index];
            int number_of_targets_in_cluster = target_end-target_start;
            int streamID = j%3;

    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/

            if (run_params->kernel == COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_Coulomb_Direct(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Coulomb_SS_Direct(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /***************** Yukawa **********************/
    /***********************************************/

            } else if (run_params->kernel == YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_Yukawa_Direct(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Yukawa_SS_Direct(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /*********** Regularized Coulomb ***************/
    /***********************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedCoulomb_Direct(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedCoulomb_SS_Direct(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /*********** Regularized Yukawa ****************/
    /***********************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedYukawa_Direct(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedYukawa_SS_Direct(number_of_targets_in_cluster, numberOfSources,
                            target_start, batchStart,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

            } else {
                printf("**ERROR** INVALID KERNEL. EXITING.\n");
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
