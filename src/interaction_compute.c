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
#include "struct_nodes.h"
#include "struct_particles.h"

#include "kernels/coulomb.h"
#include "kernels/yukawa.h"
#include "kernels/coulomb_singularity_subtraction.h"
#include "kernels/yukawa_singularity_subtraction.h"
#include "kernels/regularized-coulomb.h"
#include "kernels/regularized-yukawa.h"
#include "kernels/regularized-coulomb-singularity-subtraction.h"
#include "kernels/regularized-yukawa-singularity-subtraction.h"
#include "kernels/atan.h"



#include "interaction_compute.h"


void Interaction_PC_Compute(struct tnode_array *tree_array, struct tnode_array *batches,
                            int *tree_inter_list, int *direct_inter_list,
                            double *source_x, double *source_y, double *source_z,
                            double *source_charge, double *source_weight,
                            double *target_x, double *target_y, double *target_z, double *target_charge,
                            double *cluster_x, double *cluster_y, double *cluster_z,
                            double *cluster_charge, double *cluster_weight,
                            double *pointwisePotential, int interpolationOrder,
                            int numSources, int numTargets, int totalNumberOfInterpolationPoints,
                            int batch_approx_offset, int batch_direct_offset,
                            struct kernel *kernel, char *singularityHandling,
                            char *approximationName)
{

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

    if (strcmp(approximationName, "hermite") == 0)
        numberOfClusterCharges = 8 * totalNumberOfInterpolationPoints;

    if ((strcmp(approximationName, "hermite") == 0) && (strcmp(singularityHandling, "subtraction") == 0))
        numberOfClusterWeights = 8 * totalNumberOfInterpolationPoints;

#ifdef OPENACC_ENABLED
    #pragma acc data copyin(source_x[0:numSources], source_y[0:numSources], source_z[0:numSources], \
                        source_charge[0:numSources], source_weight[0:numSources], \
                        target_x[0:numTargets], target_y[0:numTargets], target_z[0:numTargets], \
                        target_charge[0:numTargets], \
                        cluster_x[0:totalNumberOfInterpolationPoints], cluster_y[0:totalNumberOfInterpolationPoints], \
                        cluster_z[0:totalNumberOfInterpolationPoints], \
                        cluster_charge[0:numberOfClusterCharges], cluster_weight[0:numberOfClusterWeights], \
                        tree_inter_list[0:batch_approx_offset*batch_numnodes], \
                        direct_inter_list[0:batch_direct_offset*batch_numnodes], \
                        ibegs[0:tree_numnodes], iends[0:tree_numnodes]) \
                        copy(potentialDueToApprox[0:numTargets], potentialDueToDirect[0:numTargets])
#endif
    {

    int numberOfInterpolationPoints = (interpolationOrder+1)*(interpolationOrder+1)*(interpolationOrder+1);

    for (int i = 0; i < batches->numnodes; i++) {
        int batch_ibeg = batches->ibeg[i];
        int batch_iend = batches->iend[i];
        int numberOfClusterApproximations = batches->numApprox[i];
        int numberOfDirectSums = batches->numDirect[i];

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
            if (strcmp(kernel->name, "coulomb") == 0) {

                if (strcmp(approximationName, "lagrange") == 0) {

                    if (strcmp(singularityHandling, "skipping") == 0) {
            
                        coulombApproximationLagrange(numberOfTargets,
                                numberOfInterpolationPoints, batchStart, clusterStart,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_charge,
                                kernel, potentialDueToApprox, streamID);

                    } else if (strcmp(singularityHandling, "subtraction") == 0) {

                        coulombSingularitySubtractionApproximationLagrange(numberOfTargets,
                                numberOfInterpolationPoints, batchStart, clusterStart,
                                target_x, target_y, target_z, target_charge,
                                cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                kernel, potentialDueToApprox, streamID);

                    } else {
                        printf("Invalid choice of singularityHandling. Exiting. \n");
                        exit(1);
                    }

                } else if (strcmp(approximationName, "hermite") == 0) {

                    if (strcmp(singularityHandling, "skipping") == 0) {

                        coulombApproximationHermite(numberOfTargets,
                                numberOfInterpolationPoints, batchStart,
                                clusterStart, totalNumberOfInterpolationPoints,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_charge,
                                kernel, potentialDueToApprox, streamID);

                    } else if (strcmp(singularityHandling, "subtraction") == 0) {

                        coulombSingularitySubtractionApproximationHermite(numberOfTargets,
                                numberOfInterpolationPoints, batchStart,
                                clusterStart, totalNumberOfInterpolationPoints,
                                target_x, target_y, target_z, target_charge,
                                cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                kernel, potentialDueToApprox, streamID);

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

                        yukawaApproximationLagrange(numberOfTargets,
                                numberOfInterpolationPoints, batchStart, clusterStart,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_charge,
                                kernel, potentialDueToApprox, streamID);

                    } else if (strcmp(singularityHandling, "subtraction") == 0) {
                        
                        yukawaSingularitySubtractionApproximationLagrange(numberOfTargets,
                                numberOfInterpolationPoints, batchStart, clusterStart,
                                target_x, target_y, target_z, target_charge,
                                cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                kernel, potentialDueToApprox, streamID);

                    } else {
                        printf("Invalid choice of singularityHandling. Exiting. \n");
                        exit(1);
                    }

                } else if (strcmp(approximationName, "hermite") == 0) {

                    if (strcmp(singularityHandling, "skipping") == 0) {

                        yukawaApproximationHermite(numberOfTargets,
                                numberOfInterpolationPoints, batchStart,
                                clusterStart, totalNumberOfInterpolationPoints,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_charge,
                                kernel, potentialDueToApprox, streamID);

                    } else if (strcmp(singularityHandling, "subtraction") == 0) {

                        yukawaSingularitySubtractionApproximationHermite(numberOfTargets,
                                numberOfInterpolationPoints, batchStart,
                                clusterStart, totalNumberOfInterpolationPoints,
                                target_x, target_y, target_z, target_charge,
                                cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                kernel, potentialDueToApprox, streamID);

                    } else {
                        printf("Invalid choice of singularityHandling. Exiting. \n");
                        exit(1);
                    }

                } else {
                    printf("Invalid approximationName.\n");
                    exit(1);
                }

    /***************************************/
    /********* Regularized-Coulomb *********/
    /***************************************/

            } else if (strcmp(kernel->name, "regularized-coulomb") == 0) {

                if (strcmp(approximationName, "lagrange") == 0) {

                    if (strcmp(singularityHandling, "skipping") == 0) {

                        regularizedCoulombApproximationLagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_charge,
                                    kernel, potentialDueToApprox, streamID);

                    } else if (strcmp(singularityHandling, "subtraction") == 0) {

                        regularizedCoulombSingularitySubtractionApproximationLagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z, target_charge,
                                    cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                    kernel, potentialDueToApprox, streamID);
                    }
                } else if (strcmp(approximationName, "hermite") == 0) {
                    printf("NOT SET UP FOR REGULARIZED COULOMB HERMITE.  EXITING.\n");
                    exit(-1);
                }

    /***************************************/
    /********* Regularized-Yukawa *********/
    /***************************************/

            } else if (strcmp(kernel->name, "regularized-yukawa") == 0) {

                if (strcmp(approximationName, "lagrange") == 0) {

                    if (strcmp(singularityHandling, "skipping") == 0) {

                        regularizedYukawaApproximationLagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_charge,
                                    kernel, potentialDueToApprox, streamID);

                    } else if (strcmp(singularityHandling, "subtraction") == 0) {

                        regularizedYukawaSingularitySubtractionApproximationLagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z, target_charge,
                                    cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight,
                                    kernel, potentialDueToApprox, streamID);

                    }

                } else if (strcmp(approximationName, "hermite") == 0) {
                    printf("NOT SET UP FOR REGULARIZED COULOMB HERMITE.  EXITING.\n");
                    exit(-1);
                }

    /***************************************/
    /****************** ATAN ***************/
    /***************************************/
            } else if (strcmp(kernel->name, "atan") == 0) {

                if (strcmp(approximationName, "lagrange") == 0) {

                        atanApproximationLagrange(numberOfTargets,
                                    numberOfInterpolationPoints, batchStart, clusterStart,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_charge,
                                    kernel, potentialDueToApprox, streamID);

                } else if (strcmp(approximationName, "hermite") == 0) {
                    printf("NOT SET UP FOR REGULARIZED ATAN HERMITE.  EXITING.\n");
                    exit(-1);
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
            int source_start = ibegs[node_index] - 1;
            int source_end = iends[node_index];
            int number_of_sources_in_cluster = source_end-source_start;
            int streamID = j%3;

    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/

            if (strcmp(kernel->name, "coulomb") == 0) {

                if (strcmp(singularityHandling, "skipping") == 0) {

                    coulombDirect(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_charge, source_weight,
                            kernel, potentialDueToDirect, streamID);

                } else if (strcmp(singularityHandling, "subtraction") == 0) {

                    coulombSingularitySubtractionDirect(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z, target_charge,
                            source_x, source_y, source_z, source_charge, source_weight,
                            kernel, potentialDueToDirect, streamID);

                }else {
                    printf("Invalid choice of singularityHandling. Exiting. \n");
                    exit(1);
                }

    /***********************************************/
    /***************** Yukawa **********************/
    /***********************************************/

            } else if (strcmp(kernel->name, "yukawa") == 0) {

                if (strcmp(singularityHandling, "skipping") == 0) {

                    yukawaDirect(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_charge, source_weight,
                            kernel, potentialDueToDirect, streamID);

                } else if (strcmp(singularityHandling, "subtraction") == 0) {

                    yukawaSingularitySubtractionDirect(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z, target_charge,
                            source_x, source_y, source_z, source_charge, source_weight,
                            kernel, potentialDueToDirect, streamID);

                } else {
                    printf("Invalid choice of singularityHandling. Exiting. \n");
                    exit(1);
                }

    /***************************************/
    /********* Regularized-Coulomb *********/
    /***************************************/

            } else if (strcmp(kernel->name, "regularized-coulomb") == 0) {

                if (strcmp(singularityHandling, "skipping") == 0) {

                    regularizedCoulombDirect(numberOfTargets, number_of_sources_in_cluster,
                                batchStart, source_start,
                                target_x, target_y, target_z,
                                source_x, source_y, source_z, source_charge, source_weight,
                                kernel, potentialDueToDirect, streamID);

                } else if (strcmp(singularityHandling, "subtraction") == 0) {

                    regularizedCoulombSingularitySubtractionDirect(numberOfTargets, number_of_sources_in_cluster,
                                batchStart, source_start,
                                target_x, target_y, target_z, target_charge,
                                source_x, source_y, source_z, source_charge, source_weight,
                                kernel, potentialDueToDirect, streamID);
                }


    /***************************************/
    /********* Regularized-Yukawa *********/
    /***************************************/

            } else if (strcmp(kernel->name, "regularized-yukawa") == 0) {

                if (strcmp(singularityHandling, "skipping") == 0) {

                regularizedYukawaDirect(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_charge, source_weight,
                            kernel, potentialDueToDirect, streamID);

                } else if (strcmp(singularityHandling, "subtraction") == 0) {

                    regularizedYukawaSingularitySubtractionDirect(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z, target_charge,
                            source_x, source_y, source_z, source_charge, source_weight,
                            kernel, potentialDueToDirect, streamID);

                }



            } else if (strcmp(kernel->name, "atan") == 0) {


                atanDirect(numberOfTargets, number_of_sources_in_cluster,
                            batchStart, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_charge, source_weight,
                            kernel, potentialDueToDirect, streamID);


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

    for (int k = 0; k < numTargets; k++) 
        pointwisePotential[k] += potentialDueToDirect[k];

    for (int k = 0; k < numTargets; k++) 
        pointwisePotential[k] += potentialDueToApprox[k];

    free_vector(potentialDueToDirect);
    free_vector(potentialDueToApprox);

    return;

} /* END of Interaction_PC_Compute */




void Interaction_Direct_Compute(double *source_x, double *source_y, double *source_z,
                                double *source_q, double *source_w,
                                double *target_x, double *target_y, double *target_z, double *target_q,
                                double *pointwisePotential, int numSources, int numTargets,
                                struct kernel *kernel, char *singularityHandling,
                                char *approximationName)
{


#ifdef OPENACC_ENABLED
    #pragma acc data copyin(source_x[0:numSources], source_y[0:numSources], source_z[0:numSources], \
                            source_q[0:numSources], source_w[0:numSources], \
                            target_x[0:numTargets], target_y[0:numTargets], target_z[0:numTargets], \
                            target_q[0:numTargets]), copy(pointwisePotential[0:numTargets])
#endif
    {


/**********************************************************/
/**************** COMPLETE DIRECT SUM *********************/
/**********************************************************/


    /***************************************/
    /********* Coulomb *********************/
    /***************************************/

    if (strcmp(kernel->name, "coulomb") == 0) {

        if (strcmp(singularityHandling, "skipping") == 0) {

            coulombDirect(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    kernel, pointwisePotential, 0);

        } else if (strcmp(singularityHandling, "subtraction") == 0) {

            coulombSingularitySubtractionDirect(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    kernel, pointwisePotential, 0);

        }else {
            printf("Invalid choice of singularityHandling. Exiting. \n");
            exit(1);
        }

    /***************************************/
    /********* Yukawa **********************/
    /***************************************/

    } else if (strcmp(kernel->name, "yukawa") == 0) {

        if (strcmp(singularityHandling, "skipping") == 0) {

            yukawaDirect(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    kernel, pointwisePotential, 0);

        } else if (strcmp(singularityHandling, "subtraction") == 0) {

            yukawaSingularitySubtractionDirect(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    kernel, pointwisePotential, 0);

        } else {
            printf("Invalid choice of singularityHandling. Exiting. \n");
            exit(1);
        }


    /***************************************/
    /********* Regularized-Coulomb *********/
    /***************************************/

    } else if (strcmp(kernel->name, "regularized-coulomb") == 0) {

        if (strcmp(singularityHandling, "skipping") == 0) {

            regularizedCoulombDirect(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z,
                    source_x, source_y, source_z, source_q, source_w,
                    kernel, pointwisePotential, 0);

        } else if (strcmp(singularityHandling, "subtraction") == 0) {

            regularizedCoulombSingularitySubtractionDirect(numTargets, numSources, 0, 0,
                    target_x, target_y, target_z, target_q,
                    source_x, source_y, source_z, source_q, source_w,
                    kernel, pointwisePotential, 0);

        }


    /***************************************/
    /********* Regularized-Yukawa *********/
    /***************************************/

    } else if (strcmp(kernel->name, "regularized-yukawa") == 0) {

        if (strcmp(singularityHandling, "skipping") == 0) {

            regularizedYukawaDirect(numTargets, numSources, 0, 0,
                        target_x, target_y, target_z,
                        source_x, source_y, source_z, source_q, source_w,
                        kernel, pointwisePotential, 0);

        } else if (strcmp(singularityHandling, "subtraction") == 0) {

            regularizedYukawaSingularitySubtractionDirect(numTargets, numSources, 0, 0,
                        target_x, target_y, target_z, target_q,
                        source_x, source_y, source_z, source_q, source_w,
                        kernel, pointwisePotential, 0);
        }

    } else if (strcmp(kernel->name, "atan") == 0) {

            atanDirect(numTargets, numSources, 0, 0,
                        target_x, target_y, target_z,
                        source_x, source_y, source_z, source_q, source_w,
                        kernel, pointwisePotential, 0);

    } else {
        printf("Invalid kernel->name. Exiting.\n");
        exit(1);
    }

#ifdef OPENACC_ENABLED
        #pragma acc wait
#endif
    } // end acc data region

    return;
}




void Interaction_SubtractionPotentialCorrection(double *pointwisePotential, double *target_q, int numTargets,
                                  struct kernel *kernel, char *singularityHandling)
{

    if (strcmp(singularityHandling, "subtraction") == 0) {
        if (strcmp(kernel->name, "coulomb") == 0) {
            coulombSingularitySubtractionCorrection(pointwisePotential, target_q, numTargets, kernel);

        } else if (strcmp(kernel->name, "regularized-coulomb") == 0) {
            regularizedCoulombSingularitySubtractionCorrection(pointwisePotential, target_q, numTargets, kernel);

        } else if (strcmp(kernel->name, "yukawa") == 0) {
            yukawaSingularitySubtractionCorrection(pointwisePotential, target_q, numTargets, kernel);

        } else if (strcmp(kernel->name, "regularized-yukawa") == 0) {
            regularizedYukawaSingularitySubtractionCorrection(pointwisePotential, target_q, numTargets, kernel);

        }

    }

    return;
}
