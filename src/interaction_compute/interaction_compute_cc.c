#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"
#include "../interaction_lists/struct_interaction_lists.h"

#include "../kernels/coulomb/coulomb.h"
#include "../kernels/yukawa/yukawa.h"
#include "../kernels/regularized-coulomb/regularized-coulomb.h"
#include "../kernels/regularized-yukawa/regularized-yukawa.h"
#include "../kernels/sin-over-r/sin-over-r.h"

#include "interaction_compute.h"


void InteractionCompute_CC(struct Tree *source_tree, struct Tree *target_tree,
                           struct InteractionLists *interaction_list,
                           double *source_x, double *source_y, double *source_z,
                           double *source_q, double *source_w,
                           double *target_x, double *target_y, double *target_z, double *target_q,
                           double *source_cluster_x, double *source_cluster_y, double *source_cluster_z,
                           double *source_cluster_q, double *source_cluster_w,
                           double *target_cluster_x, double *target_cluster_y, double *target_cluster_z,
                           double *target_cluster_q, double *target_cluster_w,
                           double *pointwisePotential,
                           int numSources, int numTargets, int numSourceClusterPoints, int numTargetClusterPoints,
                           struct RunParams *run_params)
{

    int **approx_inter_list = interaction_list->approx_interactions;
    int **direct_inter_list = interaction_list->direct_interactions;
    
    int *num_approx = interaction_list->num_approx;
    int *num_direct = interaction_list->num_direct;

    int source_tree_numnodes = source_tree->numnodes;
    int target_tree_numnodes = target_tree->numnodes;
    
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
    
    int *source_tree_ibeg = source_tree->ibeg;
    int *source_tree_iend = source_tree->iend;
    int *source_tree_cluster_ind = source_tree->cluster_ind;
    
    int *target_tree_ibeg = target_tree->ibeg;
    int *target_tree_iend = target_tree->iend;
    int *target_tree_cluster_ind = target_tree->cluster_ind;
    
    int numSourceClusterCharges = numSourceClusterPoints;
    int numSourceClusterWeights = numSourceClusterPoints;
    
    int numTargetClusterCharges = numTargetClusterPoints;
    int numTargetClusterWeights = numTargetClusterPoints;

    if (run_params->approximation == HERMITE) {
        numSourceClusterCharges = 8 * numSourceClusterPoints;
        numTargetClusterCharges = 8 * numTargetClusterPoints;
    }

    if ((run_params->approximation == HERMITE) && (run_params->singularity == SUBTRACTION)) {
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
                            zCT[0:numTargetClusterPoints]) \
                       copy(qCT[0:numTargetClusterCharges], \
                            wCT[0:numTargetClusterWeights], \
                            pointwisePotential[0:numTargets])
#endif
    {

    int numInterpPoints = run_params->interp_pts_per_cluster;

    for (int i = 0; i < target_tree_numnodes; i++) {
        int target_ibeg = target_tree_ibeg[i];
        int target_iend = target_tree_iend[i];
        int target_cluster_start = numInterpPoints * target_tree_cluster_ind[i];

        int numTargets = target_iend - target_ibeg + 1;
        int targetStart =  target_ibeg - 1;
        
        int numClusterApprox = num_approx[i];
        int numClusterDirect = num_direct[i];


/**********************************************************/
/************** POTENTIAL FROM APPROX *********************/
/**********************************************************/

        for (int j = 0; j < numClusterApprox; j++) {
            int source_node_index = approx_inter_list[i][j];
            int source_cluster_start = numInterpPoints * source_tree_cluster_ind[source_node_index];
            int streamID = j%3;


    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/
            if (run_params->kernel == COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {
            
                        K_Coulomb_CP_Lagrange(numInterpPoints, numInterpPoints,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CC COULOMB SS. EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** CC HERMITE CURRENTLY INOPERABLE. EXITING. \n");
                    exit(1);

                    if (run_params->singularity == SKIPPING) {

                        //K_Coulomb_CP_Hermite(numInterpPoints, numInterpPoints,
                        //    source_cluster_start, target_cluster_start,
                        //    source_cluster_x, source_cluster_y, source_cluster_z,
                        //    source_cluster_q, source_cluster_w,
                        //    target_cluster_x, target_cluster_y, target_cluster_z,
                        //    target_cluster_q,
                        //    run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CC COULOMB SS. EXITING.\n");
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

                        K_Yukawa_CP_Lagrange(numInterpPoints, numInterpPoints,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CC YUKAWA SS. EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** CC HERMITE CURRENTLY INOPERABLE. EXITING. \n");
                    exit(1);

                    if (run_params->singularity == SKIPPING) {

                        //K_Yukawa_CP_Hermite(numInterpPoints, numInterpPoints,
                        //    source_cluster_start, target_cluster_start,
                        //    source_cluster_x, source_cluster_y, source_cluster_z,
                        //    source_cluster_q, source_cluster_w,
                        //    target_cluster_x, target_cluster_y, target_cluster_z,
                        //    target_cluster_q,
                        //    run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CC YUKAWA SS. EXITING.\n");
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
    /********* Regularized Coulomb *****************/
    /***********************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_CP_Lagrange(numInterpPoints, numInterpPoints,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CC REGULARIZED COULOMB SS. EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** CC HERMITE CURRENTLY INOPERABLE. EXITING. \n");
                    exit(1);

                    if (run_params->singularity == SKIPPING) {

                        //K_RegularizedCoulomb_CP_Hermite(numInterpPoints, numInterpPoints,
                        //    source_cluster_start, target_cluster_start,
                        //    source_cluster_x, source_cluster_y, source_cluster_z,
                        //    source_cluster_q, source_cluster_w,
                        //    target_cluster_x, target_cluster_y, target_cluster_z,
                        //    target_cluster_q,
                        //    run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CC REGULARIZED COULOMB SS. EXITING.\n");
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
    /********* Regularized Yukawa ******************/
    /***********************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedYukawa_CP_Lagrange(numInterpPoints, numInterpPoints,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, streamID);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CC REGULARIZED COULOMB SS. EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** CC HERMITE CURRENTLY INOPERABLE. EXITING. \n");
                    exit(1);

                    if (run_params->singularity == SKIPPING) {

                        printf("**ERROR** NOT SET UP FOR CC REGULARIZED COULOMB HERMITE. EXITING.\n");
                        exit(1);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CC REGULARIZED COULOMB SS. EXITING.\n");
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
    /********* Sin Over R **************************/
    /***********************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_CP_Lagrange(numInterpPoints, numInterpPoints,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, streamID);
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

        for (int j = 0; j < numClusterDirect; j++) {

            int source_node_index = direct_inter_list[i][j];
            int source_ibeg = source_tree_ibeg[source_node_index];
            int source_iend = source_tree_iend[source_node_index];
            
            int numSources = source_iend - source_ibeg + 1;
            int sourceStart =  source_ibeg - 1;
            int streamID = j%3;

    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/

            if (run_params->kernel == COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_Coulomb_Direct(numTargets, numSources,
                            targetStart, sourceStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Coulomb_SS_Direct(numTargets, numSources,
                            targetStart, sourceStart,
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

                    K_Yukawa_Direct(numTargets, numSources,
                            targetStart, sourceStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Yukawa_SS_Direct(numTargets, numSources,
                            targetStart, sourceStart,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /************ Regularized Coulomb **************/
    /***********************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedCoulomb_Direct(numTargets, numSources,
                            targetStart, sourceStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedCoulomb_SS_Direct(numTargets, numSources,
                            targetStart, sourceStart,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /************ Regularized Yukawa ***************/
    /***********************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedYukawa_Direct(numTargets, numSources,
                            targetStart, sourceStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedYukawa_SS_Direct(numTargets, numSources,
                            targetStart, sourceStart,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /************ Sin Over R ***********************/
    /***********************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                if (run_params->singularity == SKIPPING) {

                    K_SinOverR_Direct(numTargets, numSources,
                            targetStart, sourceStart,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, pointwisePotential, streamID);
                }

            } else {
                printf("**ERROR** INVALID KERNEL. EXITING.\n");
                exit(1);
            }

        } // end loop over number of direct interactions

    } // end loop over target nodes

#ifdef OPENACC_ENABLED
        #pragma acc wait
#endif
    } // end acc data region

    return;

} /* END of function cc_treecode */
