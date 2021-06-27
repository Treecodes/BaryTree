#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"
#include "../clusters/struct_clusters.h"
#include "../run_params/struct_run_params.h"
#include "../interaction_lists/struct_interaction_lists.h"

#include "../kernels/coulomb/coulomb.h"
#include "../kernels/yukawa/yukawa.h"
#include "../kernels/regularized-coulomb/regularized-coulomb.h"
#include "../kernels/regularized-yukawa/regularized-yukawa.h"
#include "../kernels/sin-over-r/sin-over-r.h"
#include "../kernels/rbs-u/rbs-u.h"
#include "../kernels/rbs-v/rbs-v.h"
#include "../kernels/user_kernel/user_kernel.h"

#include "interaction_compute.h"


void InteractionCompute_CC(double *potential, struct Tree *source_tree, struct Tree *target_tree,
                           struct InteractionLists *interaction_list,
                           struct Particles *sources, struct Particles *targets,
                           struct Clusters *source_clusters, struct Clusters *target_clusters,
                           struct RunParams *run_params)
{
    int interp_pts_per_cluster = run_params->interp_pts_per_cluster;

    int **direct_inter_list = interaction_list->pp_interactions;
    int **approx_inter_list = interaction_list->cc_interactions;
    int **source_approx_inter_list = interaction_list->pc_interactions;
    int **target_approx_inter_list = interaction_list->cp_interactions;

    int *num_direct = interaction_list->num_pp;
    int *num_approx = interaction_list->num_cc;
    int *num_source_approx = interaction_list->num_pc;
    int *num_target_approx = interaction_list->num_cp;

    int source_tree_numnodes = source_tree->numnodes;
    int target_tree_numnodes = target_tree->numnodes;
    
    int num_sources  = sources->num;
    double *source_x = sources->x;
    double *source_y = sources->y;
    double *source_z = sources->z;
    double *source_q = sources->q;
    double *source_w = sources->w;

    int num_targets  = targets->num;
    double *target_x = targets->x;
    double *target_y = targets->y;
    double *target_z = targets->z;
    double *target_q = targets->q;

    int num_source_cluster_points  = source_clusters->num;
    int num_source_cluster_charges = source_clusters->num_charges;
    int num_source_cluster_weights = source_clusters->num_weights;
    double *source_cluster_x = source_clusters->x;
    double *source_cluster_y = source_clusters->y;
    double *source_cluster_z = source_clusters->z;
    double *source_cluster_q = source_clusters->q;
    double *source_cluster_w = source_clusters->w;
    
    int num_target_cluster_points  = target_clusters->num;
    int num_target_cluster_charges = target_clusters->num_charges;
    int num_target_cluster_weights = target_clusters->num_weights;
    double *target_cluster_x = target_clusters->x;
    double *target_cluster_y = target_clusters->y;
    double *target_cluster_z = target_clusters->z;
    double *target_cluster_q = target_clusters->q;
    double *target_cluster_w = target_clusters->w;
    
    int *source_tree_ibeg = source_tree->ibeg;
    int *source_tree_iend = source_tree->iend;
    int *source_tree_cluster_ind = source_tree->cluster_ind;
    
    int *target_tree_ibeg = target_tree->ibeg;
    int *target_tree_iend = target_tree->iend;
    int *target_tree_cluster_ind = target_tree->cluster_ind;
    
    // Additionally, not setup for Hermite either at the moment.
        
    for (int i = 0; i < target_tree_numnodes; i++) {
        int target_ibeg = target_tree_ibeg[i];
        int target_iend = target_tree_iend[i];
        int target_cluster_start = interp_pts_per_cluster * target_tree_cluster_ind[i];

        int num_targets_in_cluster = target_iend - target_ibeg + 1;
        int target_start =  target_ibeg - 1;
        
        int num_approx_in_cluster = num_approx[i];
        int num_direct_in_cluster = num_direct[i];
        
        int num_source_approx_in_cluster = num_source_approx[i];
        int num_target_approx_in_cluster = num_target_approx[i];


/* * ********************************************************/
/* * ************ POTENTIAL FROM APPROX *********************/
/* * ********************************************************/

        for (int j = 0; j < num_approx_in_cluster; j++) {
            int source_node_index = approx_inter_list[i][j];
            int source_cluster_start = interp_pts_per_cluster * source_tree_cluster_ind[source_node_index];
            int stream_id = j%3;


    /* * *********************************************/
    /* * *************** Coulomb *********************/
    /* * *********************************************/
            if (run_params->kernel == COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {
            
                        K_Coulomb_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Coulomb_SS_CC_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q, target_cluster_w,
                            run_params, stream_id);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** CC HERMITE CURRENTLY INOPERABLE. EXITING. \n");
                    exit(1);

                    if (run_params->singularity == SKIPPING) {

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

    /* * *********************************************/
    /* * *************** Yukawa **********************/
    /* * *********************************************/

            } else if (run_params->kernel == YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Yukawa_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Yukawa_SS_CC_Lagrange(interp_pts_per_cluster,
                            interp_pts_per_cluster, source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q, target_cluster_w,
                            run_params, stream_id);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** CC HERMITE CURRENTLY INOPERABLE. EXITING. \n");
                    exit(1);

                    if (run_params->singularity == SKIPPING) {

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

    /* * *********************************************/
    /* * ******* Regularized Coulomb *****************/
    /* * *********************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);

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

    /* * *********************************************/
    /* * ******* Regularized Yukawa ******************/
    /* * *********************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedYukawa_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);

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

    /* * *********************************************/
    /* * ******* Sin Over R **************************/
    /* * *********************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);
                    }
                }

    /* * *********************************************/
    /* * ******* RBS_U *******************************/
    /* * *********************************************/

            } else if (run_params->kernel == RBS_U) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RBSu_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);
                    }
                }

    /* * *********************************************/
    /* * ******* RBS_V *******************************/
    /* * *********************************************/

            } else if (run_params->kernel == RBS_V) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RBSv_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                            source_cluster_start, target_cluster_start,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);
                    }
                }

    /* * *********************************************/
    /* * ******* USER DEFINED KERNEL *****************/
    /* * *********************************************/

            } else if (run_params->kernel == USER) {

                if (run_params->approximation == LAGRANGE) {

                    K_User_Kernel_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                        source_cluster_start, target_cluster_start,
                        source_cluster_x, source_cluster_y, source_cluster_z,
                        source_cluster_q,
                        target_cluster_x, target_cluster_y, target_cluster_z,
                        target_cluster_q,
                        run_params, stream_id);

                }

            } else {
                printf("**ERROR** INVALID KERNEL. EXITING.\n");
                exit(1);
            }

        } // end loop over cluster approximations
        
        
        
/* * ********************************************************/
/* * ************ POTENTIAL FROM SOURCE APPROX (PC) *********/
/* * ********************************************************/

        for (int j = 0; j < num_source_approx_in_cluster; j++) {
            int source_node_index = source_approx_inter_list[i][j];
            int source_cluster_start = interp_pts_per_cluster * source_tree_cluster_ind[source_node_index];
            int stream_id = j%3;

    /* * *********************************************/
    /* * *************** Coulomb *********************/
    /* * *********************************************/
            if (run_params->kernel == COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {
            
                        K_Coulomb_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                            target_start, source_cluster_start,
                            target_x, target_y, target_z,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            run_params, potential, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Coulomb_SS_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                            target_start, source_cluster_start,
                            target_x, target_y, target_z, target_q,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            run_params, potential, stream_id);


                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** CC HERMITE CURRENTLY INOPERABLE. EXITING. \n");
                    exit(1);

                    if (run_params->singularity == SKIPPING) {

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

    /* * *********************************************/
    /* * *************** Yukawa **********************/
    /* * *********************************************/

            } else if (run_params->kernel == YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Yukawa_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                            target_start, source_cluster_start,
                            target_x, target_y, target_z,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            run_params, potential, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Yukawa_SS_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                            target_start, source_cluster_start,
                            target_x, target_y, target_z, target_q,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q, source_cluster_w,
                            run_params, potential, stream_id);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** CC HERMITE CURRENTLY INOPERABLE. EXITING. \n");
                    exit(1);

                    if (run_params->singularity == SKIPPING) {

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

    /* * *********************************************/
    /* * ******* Regularized Coulomb *****************/
    /* * *********************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                            target_start, source_cluster_start,
                            target_x, target_y, target_z,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            run_params, potential, stream_id);

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

    /* * *********************************************/
    /* * ******* Regularized Yukawa ******************/
    /* * *********************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedYukawa_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                            target_start, source_cluster_start,
                            target_x, target_y, target_z,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            run_params, potential, stream_id);

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

    /* * *********************************************/
    /* * ******* Sin Over R **************************/
    /* * *********************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                            target_start, source_cluster_start,
                            target_x, target_y, target_z,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            run_params, potential, stream_id);
                    }
                }

    /* * *********************************************/
    /* * ******* RBS U *******************************/
    /* * *********************************************/

            } else if (run_params->kernel == RBS_U) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RBSu_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                            target_start, source_cluster_start,
                            target_x, target_y, target_z,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            run_params, potential, stream_id);
                    }
                }

    /* * *********************************************/
    /* * ******* RBS V *******************************/
    /* * *********************************************/

            } else if (run_params->kernel == RBS_V) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RBSv_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                            target_start, source_cluster_start,
                            target_x, target_y, target_z,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            run_params, potential, stream_id);
                    }
                }

    /* * *********************************************/
    /* * ******* USER DEFINED KERNEL *****************/
    /* * *********************************************/

            } else if (run_params->kernel == USER) {

                if (run_params->approximation == LAGRANGE) {

                        K_User_Kernel_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                            target_start, source_cluster_start,
                            target_x, target_y, target_z,
                            source_cluster_x, source_cluster_y, source_cluster_z,
                            source_cluster_q,
                            run_params, potential, stream_id);

                }

            } else {
                printf("[Interaction_Compute_CC] **ERROR** INVALID KERNEL. EXITING.\n");
                exit(1);
            }

        } // end loop over cluster approximations



/* * ********************************************************/
/* * ************ POTENTIAL FROM TARGET APPROX (CP) *********/
/* * ********************************************************/

        for (int j = 0; j < num_target_approx_in_cluster; j++) {
        
            int source_node_index = target_approx_inter_list[i][j];
            int source_ibeg = source_tree_ibeg[source_node_index];
            int source_iend = source_tree_iend[source_node_index];
            
            int num_sources_in_cluster = source_iend - source_ibeg + 1;
            int source_start =  source_ibeg - 1;
            int stream_id = j%3;

    /* * *********************************************/
    /* * *************** Coulomb *********************/
    /* * *********************************************/
            if (run_params->kernel == COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {
            
                        K_Coulomb_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                            source_start, target_cluster_start,
                            source_x, source_y, source_z, source_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Coulomb_SS_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                            source_start, target_cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q, target_cluster_w,
                            run_params, stream_id);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** CC HERMITE CURRENTLY INOPERABLE. EXITING. \n");
                    exit(1);

                    if (run_params->singularity == SKIPPING) {

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

    /* * *********************************************/
    /* * *************** Yukawa **********************/
    /* * *********************************************/

            } else if (run_params->kernel == YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Yukawa_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                            source_start, target_cluster_start,
                            source_x, source_y, source_z, source_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Yukawa_SS_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                            source_start, target_cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q, target_cluster_w,
                            run_params, stream_id);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** CC HERMITE CURRENTLY INOPERABLE. EXITING. \n");
                    exit(1);

                    if (run_params->singularity == SKIPPING) {

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

    /* * *********************************************/
    /* * ******* Regularized Coulomb *****************/
    /* * *********************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                            source_start, target_cluster_start,
                            source_x, source_y, source_z, source_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);

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

    /* * *********************************************/
    /* * ******* Regularized Yukawa ******************/
    /* * *********************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedYukawa_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                            source_start, target_cluster_start,
                            source_x, source_y, source_z, source_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);

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

    /* * *********************************************/
    /* * ******* Sin Over R **************************/
    /* * *********************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                            source_start, target_cluster_start,
                            source_x, source_y, source_z, source_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);
                    }
                }

    /* * *********************************************/
    /* * ******* RBS U *******************************/
    /* * *********************************************/

            } else if (run_params->kernel == RBS_U) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RBSu_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                            source_start, target_cluster_start,
                            source_x, source_y, source_z, source_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);
                    }
                }

    /* * *********************************************/
    /* * ******* RBS V *******************************/
    /* * *********************************************/

            } else if (run_params->kernel == RBS_V) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RBSv_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                            source_start, target_cluster_start,
                            source_x, source_y, source_z, source_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);
                    }
                }

    /* * *********************************************/
    /* * ******* USER DEFINED KERNEL *****************/
    /* * *********************************************/

            } else if (run_params->kernel == USER) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_User_Kernel_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                            source_start, target_cluster_start,
                            source_x, source_y, source_z, source_q,
                            target_cluster_x, target_cluster_y, target_cluster_z,
                            target_cluster_q,
                            run_params, stream_id);
                    }
                }

            } else {
                printf("**ERROR** INVALID KERNEL. EXITING.\n");
                exit(1);
            }

        } // end loop over cluster approximations



/* * ********************************************************/
/* * ************ POTENTIAL FROM DIRECT *********************/
/* * ********************************************************/

        for (int j = 0; j < num_direct_in_cluster; j++) {

            int source_node_index = direct_inter_list[i][j];
            int source_ibeg = source_tree_ibeg[source_node_index];
            int source_iend = source_tree_iend[source_node_index];
            
            int num_sources_in_cluster = source_iend - source_ibeg + 1;
            int source_start =  source_ibeg - 1;
            int stream_id = j%3;

    /* * *********************************************/
    /* * *************** Coulomb *********************/
    /* * *********************************************/

            if (run_params->kernel == COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_Coulomb_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Coulomb_SS_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /* * *********************************************/
    /* * *************** Yukawa **********************/
    /* * *********************************************/

            } else if (run_params->kernel == YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_Yukawa_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Yukawa_SS_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /* * *********************************************/
    /* * ********** Regularized Coulomb **************/
    /* * *********************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedCoulomb_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedCoulomb_SS_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /* * *********************************************/
    /* * ********** Regularized Yukawa ***************/
    /* * *********************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedYukawa_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedYukawa_SS_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /* * *********************************************/
    /* * ********** Sin Over R ***********************/
    /* * *********************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                if (run_params->singularity == SKIPPING) {

                    K_SinOverR_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);
                }

    /* * *********************************************/
    /* * ********** RBS U ****************************/
    /* * *********************************************/

            } else if (run_params->kernel == RBS_U) {

                if (run_params->singularity == SKIPPING) {

                    K_RBSu_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);
                }

    /* * *********************************************/
    /* * ********** RBS V ****************************/
    /* * *********************************************/

            } else if (run_params->kernel == RBS_V) {

                if (run_params->singularity == SKIPPING) {

                    K_RBSv_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);
                }

    /* * *********************************************/
    /* * ********** USER DEFINED KERNEL **************/
    /* * *********************************************/

            } else if (run_params->kernel == USER) {

                if (run_params->singularity == SKIPPING) {

                    K_User_Kernel_PP(num_targets_in_cluster, num_sources_in_cluster,
                            target_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);
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

    return;

} /* END of function cc_treecode */
