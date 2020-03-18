/*
 *Procedures for Particle-Cluster Treecode
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
#include "../interaction_lists/struct_interaction_lists.h"

#include "../kernels/coulomb/coulomb.h"
#include "../kernels/yukawa/yukawa.h"
#include "../kernels/regularized-coulomb/regularized-coulomb.h"
#include "../kernels/regularized-yukawa/regularized-yukawa.h"
#include "../kernels/atan/atan.h"
#include "../kernels/sin-over-r/sin-over-r.h"


#include "interaction_compute.h"


void InteractionCompute_PC(struct tnode_array *tree_array, struct tnode_array *batches,
                           struct InteractionLists *interaction_list,
                           double *source_x, double *source_y, double *source_z,
                           double *source_charge, double *source_weight,
                           double *target_x, double *target_y, double *target_z, double *target_charge,
                           double *cluster_x, double *cluster_y, double *cluster_z,
                           double *cluster_charge, double *cluster_weight,
                           double *pointwisePotential,
                           int numSources, int numTargets, int totalNumberOfInterpolationPoints,
                           struct RunParams *run_params)
{

    int **approx_inter_list = interaction_list->approx_interactions;
    int **direct_inter_list = interaction_list->direct_interactions;
    
    int *num_approx = interaction_list->num_approx;
    int *num_direct = interaction_list->num_direct;

    int tree_numnodes = tree_array->numnodes;
    int batch_numnodes = batches->numnodes;

    double *potentialDueToDirect, *potentialDueToApprox;
    make_vector(potentialDueToDirect, numTargets);
    make_vector(potentialDueToApprox, numTargets);

    memset(potentialDueToApprox, 0, numTargets * sizeof(double));
    memset(potentialDueToDirect, 0, numTargets * sizeof(double));

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
    #pragma acc data copyin(source_x[0:numSources], source_y[0:numSources], source_z[0:numSources], \
                        source_charge[0:numSources], source_weight[0:numSources], \
                        target_x[0:numTargets], target_y[0:numTargets], target_z[0:numTargets], \
                        target_charge[0:numTargets], \
                        cluster_x[0:totalNumberOfInterpolationPoints], cluster_y[0:totalNumberOfInterpolationPoints], \
                        cluster_z[0:totalNumberOfInterpolationPoints], \
                        cluster_charge[0:numberOfClusterCharges], cluster_weight[0:numberOfClusterWeights]) \
                        copy(potentialDueToApprox[0:numTargets], potentialDueToDirect[0:numTargets])
#endif
    {

    int numberOfInterpolationPoints = run_params->interp_pts_per_cluster;

    for (int i = 0; i < batches->numnodes; i++) {
        int batch_ibeg = batches->ibeg[i];
        int batch_iend = batches->iend[i];
        int numberOfClusterApproximations = num_approx[i];
        int numberOfDirectSums = num_direct[i];

        int numberOfTargets = batch_iend - batch_ibeg + 1;
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
            
                        K_Coulomb_PC_Lagrange(numberOfTargets,
                                numberOfInterpolationPoints, batchStart, clusterStart,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_charge,
                                run_params, potentialDueToApprox, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Coulomb_SS_PC_Lagrange(numberOfTargets,
                                numberOfInterpolationPoints, batchStart, clusterStart,
                                target_x, target_y, target_z, target_charge,
                                cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                run_params, potentialDueToApprox, streamID);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Coulomb_PC_Hermite(numberOfTargets,
                                numberOfInterpolationPoints, batchStart,
                                clusterStart, totalNumberOfInterpolationPoints,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_charge,
                                run_params, potentialDueToApprox, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Coulomb_SS_PC_Hermite(numberOfTargets,
                                numberOfInterpolationPoints, batchStart,
                                clusterStart, totalNumberOfInterpolationPoints,
                                target_x, target_y, target_z, target_charge,
                                cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                run_params, potentialDueToApprox, streamID);

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

                        K_Yukawa_PC_Lagrange(numberOfTargets,
                                numberOfInterpolationPoints, batchStart, clusterStart,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_charge,
                                run_params, potentialDueToApprox, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {
                        
                        K_Yukawa_SS_PC_Lagrange(numberOfTargets,
                                numberOfInterpolationPoints, batchStart, clusterStart,
                                target_x, target_y, target_z, target_charge,
                                cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                run_params, potentialDueToApprox, streamID);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Yukawa_PC_Hermite(numberOfTargets,
                                numberOfInterpolationPoints, batchStart,
                                clusterStart, totalNumberOfInterpolationPoints,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_charge,
                                run_params, potentialDueToApprox, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Yukawa_SS_PC_Hermite(numberOfTargets,
                                numberOfInterpolationPoints, batchStart,
                                clusterStart, totalNumberOfInterpolationPoints,
                                target_x, target_y, target_z, target_charge,
                                cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                run_params, potentialDueToApprox, streamID);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else {
                    printf("**ERROR** INVALID CHOICE OF APPROXIMATION. EXITING. \n");
                    exit(1);
                }

    /***************************************/
    /********* Regularized-Coulomb *********/
    /***************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_PC_Lagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_charge,
                                    run_params, potentialDueToApprox, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_RegularizedCoulomb_SS_PC_Lagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z, target_charge,
                                    cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                    run_params, potentialDueToApprox, streamID);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_PC_Hermite(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    totalNumberOfInterpolationPoints, target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_charge,
                                    run_params, potentialDueToApprox, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_RegularizedCoulomb_SS_PC_Hermite(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    totalNumberOfInterpolationPoints, target_x, target_y, target_z, target_charge,
                                    cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                    run_params, potentialDueToApprox, streamID);
                    }
                }

    /***************************************/
    /********* Regularized-Yukawa *********/
    /***************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedYukawa_PC_Lagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_charge,
                                    run_params, potentialDueToApprox, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_RegularizedYukawa_SS_PC_Lagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z, target_charge,
                                    cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                    run_params, potentialDueToApprox, streamID);

                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        printf("**ERROR** NOT SET UP FOR PC REGULARIZED YUKAWA HERMITE.  EXITING.\n");
                        exit(1);

                        /*
                        K_RegularizedYukawa_PC_Hermite(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    totalNumberOfInterpolationPoints, target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_charge,
                                    run_params, potentialDueToApprox, streamID);
                        */

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR "
                               "PC REGULARIZED YUKAWA SINGULARITY SUBTRACTION HERMITE.  EXITING.\n");
                        exit(1);

                        /*
                        K_RegularizedYukawa_SS_PC_Hermite(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    totalNumberOfInterpolationPoints, target_x, target_y, target_z, target_charge,
                                    cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                    run_params, potentialDueToApprox, streamID);
                        */
                    }
                }


    /***************************************/
    /****************** Atan ***************/
    /***************************************/

            } else if (run_params->kernel == ATAN) {

                if (run_params->approximation == LAGRANGE) {

                        K_Atan_PC_Lagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_charge,
                                    run_params, potentialDueToApprox, streamID);

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** NOT SET UP FOR PC ATAN HERMITE.  EXITING.\n");
                    exit(1);
                }


    /***************************************/
    /************* Sin Over R **************/
    /***************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_PC_Lagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_charge,
                                    run_params, potentialDueToApprox, streamID);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_PC_Hermite(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    totalNumberOfInterpolationPoints, target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_charge,
                                    run_params, potentialDueToApprox, streamID);
                    }
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
            int source_start = ibegs[node_index] - 1;
            int source_end = iends[node_index];
            int number_of_sources_in_cluster = source_end-source_start;
            int streamID = j%3;

    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/

            if (run_params->kernel == COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_Coulomb_Direct(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_charge, source_weight,
                            run_params, potentialDueToDirect, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Coulomb_SS_Direct(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z, target_charge,
                            source_x, source_y, source_z, source_charge, source_weight,
                            run_params, potentialDueToDirect, streamID);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /***************** Yukawa **********************/
    /***********************************************/

            } else if (run_params->kernel == YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_Yukawa_Direct(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_charge, source_weight,
                            run_params, potentialDueToDirect, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Yukawa_SS_Direct(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z, target_charge,
                            source_x, source_y, source_z, source_charge, source_weight,
                            run_params, potentialDueToDirect, streamID);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***************************************/
    /********* Regularized-Coulomb *********/
    /***************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedCoulomb_Direct(numberOfTargets, number_of_sources_in_cluster,
                                batchStart, source_start,
                                target_x, target_y, target_z,
                                source_x, source_y, source_z, source_charge, source_weight,
                                run_params, potentialDueToDirect, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedCoulomb_SS_Direct(numberOfTargets, number_of_sources_in_cluster,
                                batchStart, source_start,
                                target_x, target_y, target_z, target_charge,
                                source_x, source_y, source_z, source_charge, source_weight,
                                run_params, potentialDueToDirect, streamID);
                }


    /***************************************/
    /********* Regularized-Yukawa **********/
    /***************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedYukawa_Direct(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_charge, source_weight,
                            run_params, potentialDueToDirect, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedYukawa_SS_Direct(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z, target_charge,
                            source_x, source_y, source_z, source_charge, source_weight,
                            run_params, potentialDueToDirect, streamID);

                }


    /***************************************/
    /********* Atan ************************/
    /***************************************/

            } else if (run_params->kernel == ATAN) {

                K_Atan_Direct(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_charge, source_weight,
                            run_params, potentialDueToDirect, streamID);


    /***************************************/
    /********* Sin Over R ******************/
    /***************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                K_SinOverR_Direct(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_charge, source_weight,
                            run_params, potentialDueToDirect, streamID);

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

    for (int k = 0; k < numTargets; k++) 
        pointwisePotential[k] += potentialDueToDirect[k];

    for (int k = 0; k < numTargets; k++) 
        pointwisePotential[k] += potentialDueToApprox[k];

    free_vector(potentialDueToDirect);
    free_vector(potentialDueToApprox);

    return;

} /* END of Interaction_PC_Compute */
