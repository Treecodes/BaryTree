#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "../utilities/array.h"

#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"
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


void InteractionCompute_CP(double *potential, struct Tree *tree, struct Tree *batches,
                           struct InteractionLists *interaction_list,
                           struct Particles *sources, struct Particles *targets,
                           struct Clusters *clusters, struct RunParams *run_params)
{
    int interp_pts_per_cluster = run_params->interp_pts_per_cluster;

    int num_sources   = sources->num;
    double *source_x  = sources->x;
    double *source_y  = sources->y;
    double *source_z  = sources->z;
    double *source_q  = sources->q;
    double *source_w  = sources->w;

    int num_targets   = targets->num;
    double *target_x  = targets->x;
    double *target_y  = targets->y;
    double *target_z  = targets->z;
    double *target_q  = targets->q;

    int total_num_interp_pts     = clusters->num;
    int total_num_interp_charges = clusters->num_charges;
    int total_num_interp_weights = clusters->num_weights;
    double *cluster_x = clusters->x;
    double *cluster_y = clusters->y;
    double *cluster_z = clusters->z;
    double *cluster_q = clusters->q;
    double *cluster_w = clusters->w;

    int **approx_inter_list = interaction_list->cp_interactions;
    int **direct_inter_list = interaction_list->pp_interactions;
    
    int *num_approx = interaction_list->num_cp;
    int *num_direct = interaction_list->num_pp;
    
    int tree_numnodes = tree->numnodes;
    int batch_numnodes = batches->numnodes;
    
    int *tree_ibeg = tree->ibeg;
    int *tree_iend = tree->iend;
    int *cluster_ind = tree->cluster_ind;


    for (int i = 0; i < batches->numnodes; i++) {
    
        int batch_ibeg = batches->ibeg[i];
        int batch_iend = batches->iend[i];
        
        int num_approx_in_batch = num_approx[i];
        int num_direct_in_batch = num_direct[i];

        int num_sources_in_batch = batch_iend - batch_ibeg + 1;
        int batch_start =  batch_ibeg - 1;


/**********************************************************/
/************** POTENTIAL FROM APPROX *********************/
/**********************************************************/

        for (int j = 0; j < num_approx_in_batch; j++) {

            int node_index = approx_inter_list[i][j];

            int cluster_start = interp_pts_per_cluster*cluster_ind[node_index];
            int stream_id = j%3;


    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/
            if (run_params->kernel == COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {
            
                        K_Coulomb_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Coulomb_SS_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q, cluster_w,
                            run_params, stream_id);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Coulomb_CP_Hermite(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);

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

                        K_Yukawa_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Yukawa_SS_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q, cluster_w,
                            run_params, stream_id);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Yukawa_CP_Hermite(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);

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

                        K_RegularizedCoulomb_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP REGULARIZED COULOMB SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_CP_Hermite(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);

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

                        K_RegularizedYukawa_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);

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


    /***************************************/
    /********* Sin Over R ******************/
    /***************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_CP_Hermite(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);
                    }
                }

    /***************************************/
    /********* RBS U ***********************/
    /***************************************/

            } else if (run_params->kernel == RBS_U) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RBSu_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);
                    }
                }

    /***************************************/
    /********* RBS V ***********************/
    /***************************************/

            } else if (run_params->kernel == RBS_V) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RBSv_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);
                    }
                }

    /***************************************/
    /****** USER DEFINED KERNEL ************/
    /***************************************/

            } else if (run_params->kernel == USER) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_User_Kernel_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);
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

        for (int j = 0; j < num_direct_in_batch; j++) {

            int node_index = direct_inter_list[i][j];

            int target_start = tree_ibeg[node_index] - 1;
            int target_end = tree_iend[node_index];
            int num_targets_in_cluster = target_end-target_start;
            int stream_id = j%3;

    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/

            if (run_params->kernel == COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_Coulomb_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Coulomb_SS_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /***************** Yukawa **********************/
    /***********************************************/

            } else if (run_params->kernel == YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_Yukawa_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Yukawa_SS_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /*********** Regularized Coulomb ***************/
    /***********************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedCoulomb_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedCoulomb_SS_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***********************************************/
    /*********** Regularized Yukawa ****************/
    /***********************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedYukawa_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedYukawa_SS_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }


    /***********************************************/
    /*********** Sin Over R ************************/
    /***********************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                if (run_params->singularity == SKIPPING) {

                    K_SinOverR_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);
                }

    /***********************************************/
    /*********** RBS U *****************************/
    /***********************************************/

            } else if (run_params->kernel == RBS_U) {

                if (run_params->singularity == SKIPPING) {

                    K_RBSu_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);
                }

    /***********************************************/
    /*********** RBS V *****************************/
    /***********************************************/

            } else if (run_params->kernel == RBS_V) {

                if (run_params->singularity == SKIPPING) {

                    K_RBSv_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);
                }

    /***********************************************/
    /*********** User Defined Kernel ***************/
    /***********************************************/

            } else if (run_params->kernel == USER) {

                if (run_params->singularity == SKIPPING) {

                    K_User_Kernel_PP(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);
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

    return;

} /* END of function pc_treecode */
