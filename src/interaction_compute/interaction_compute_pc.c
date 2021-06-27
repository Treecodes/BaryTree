#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "../utilities/array.h"

#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"
#include "../clusters/struct_clusters.h"
#include "../run_params/struct_run_params.h"
#include "../interaction_lists/struct_interaction_lists.h"

#include "../kernels/coulomb/coulomb.h"
#include "../kernels/yukawa/yukawa.h"
#include "../kernels/regularized-coulomb/regularized-coulomb.h"
#include "../kernels/regularized-yukawa/regularized-yukawa.h"
#include "../kernels/atan/atan.h"
#include "../kernels/sin-over-r/sin-over-r.h"
#include "../kernels/mq/mq.h"
#include "../kernels/rbs-u/rbs-u.h"
#include "../kernels/rbs-v/rbs-v.h"
#include "../kernels/user_kernel/user_kernel.h"

#include "interaction_compute.h"


void InteractionCompute_PC(double *potential, struct Tree *tree, struct Tree *batches,
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

    int **approx_inter_list = interaction_list->pc_interactions;
    int **direct_inter_list = interaction_list->pp_interactions;
    
    int *num_approx = interaction_list->num_pc;
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

        int num_targets_in_batch = batch_iend - batch_ibeg + 1;
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
            
                        K_Coulomb_PC_Lagrange(num_targets_in_batch,
                                interp_pts_per_cluster, batch_start, cluster_start,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_q,
                                run_params, potential, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Coulomb_SS_PC_Lagrange(num_targets_in_batch,
                                interp_pts_per_cluster, batch_start, cluster_start,
                                target_x, target_y, target_z, target_q,
                                cluster_x, cluster_y, cluster_z, cluster_q, cluster_w,
                                run_params, potential, stream_id);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Coulomb_PC_Hermite(num_targets_in_batch,
                                interp_pts_per_cluster, batch_start,
                                cluster_start, total_num_interp_pts,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_q,
                                run_params, potential, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Coulomb_SS_PC_Hermite(num_targets_in_batch,
                                interp_pts_per_cluster, batch_start,
                                cluster_start, total_num_interp_pts,
                                target_x, target_y, target_z, target_q,
                                cluster_x, cluster_y, cluster_z, cluster_q, cluster_w,
                                run_params, potential, stream_id);

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

                        K_Yukawa_PC_Lagrange(num_targets_in_batch,
                                interp_pts_per_cluster, batch_start, cluster_start,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_q,
                                run_params, potential, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {
                        
                        K_Yukawa_SS_PC_Lagrange(num_targets_in_batch,
                                interp_pts_per_cluster, batch_start, cluster_start,
                                target_x, target_y, target_z, target_q,
                                cluster_x, cluster_y, cluster_z, cluster_q, cluster_w,
                                run_params, potential, stream_id);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Yukawa_PC_Hermite(num_targets_in_batch,
                                interp_pts_per_cluster, batch_start,
                                cluster_start, total_num_interp_pts,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_q,
                                run_params, potential, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_Yukawa_SS_PC_Hermite(num_targets_in_batch,
                                interp_pts_per_cluster, batch_start,
                                cluster_start, total_num_interp_pts,
                                target_x, target_y, target_z, target_q,
                                cluster_x, cluster_y, cluster_z, cluster_q, cluster_w,
                                run_params, potential, stream_id);

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

                        K_RegularizedCoulomb_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_RegularizedCoulomb_SS_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z, target_q,
                                    cluster_x, cluster_y, cluster_z, cluster_q, cluster_w,
                                    run_params, potential, stream_id);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_PC_Hermite(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    total_num_interp_pts, target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_RegularizedCoulomb_SS_PC_Hermite(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    total_num_interp_pts, target_x, target_y, target_z, target_q,
                                    cluster_x, cluster_y, cluster_z, cluster_q, cluster_w,
                                    run_params, potential, stream_id);
                    }
                }

    /***************************************/
    /********* Regularized-Yukawa *********/
    /***************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedYukawa_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        K_RegularizedYukawa_SS_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z, target_q,
                                    cluster_x, cluster_y, cluster_z, cluster_q, cluster_w,
                                    run_params, potential, stream_id);

                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        printf("**ERROR** NOT SET UP FOR PC REGULARIZED YUKAWA HERMITE.  EXITING.\n");
                        exit(1);

                        /*
                        K_RegularizedYukawa_PC_Hermite(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    total_num_interp_pts, target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);
                        */

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR "
                               "PC REGULARIZED YUKAWA SINGULARITY SUBTRACTION HERMITE.  EXITING.\n");
                        exit(1);

                        /*
                        K_RegularizedYukawa_SS_PC_Hermite(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    total_num_interp_pts, target_x, target_y, target_z, target_q,
                                    cluster_x, cluster_y, cluster_z, cluster_q, cluster_w,
                                    run_params, potential, stream_id);
                        */
                    }
                }


    /***************************************/
    /****************** Atan ***************/
    /***************************************/

            } else if (run_params->kernel == ATAN) {

                if (run_params->approximation == LAGRANGE) {

                        K_Atan_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);

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

                        K_SinOverR_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_PC_Hermite(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    total_num_interp_pts, target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);
                    }
                }


    /***************************************/
    /******************* MQ ****************/
    /***************************************/

            } else if (run_params->kernel == MQ) {

                if (run_params->approximation == LAGRANGE) {

                        K_MQ_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);

                } else if (run_params->approximation == HERMITE) {

                    printf("**ERROR** NOT SET UP FOR PC MQ HERMITE.  EXITING.\n");
                    exit(1);
                }

    /***************************************/
    /************* RBS U *******************/
    /***************************************/

            } else if (run_params->kernel == RBS_U) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RBSu_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);
                    }
                }

    /***************************************/
    /************* RBS V *******************/
    /***************************************/

            } else if (run_params->kernel == RBS_V) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RBSv_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);
                    }
                }

    /***************************************/
    /********* USER DEFINED KERNEL *********/
    /***************************************/

            } else if (run_params->kernel == USER) {

                if (run_params->approximation == LAGRANGE) {

                        K_User_Kernel_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential, stream_id);

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
            int source_start = tree_ibeg[node_index] - 1;
            int source_end = tree_iend[node_index];
            int num_sources_in_cluster = source_end-source_start;
            int stream_id = j%3;

    /***********************************************/
    /***************** Coulomb *********************/
    /***********************************************/

            if (run_params->kernel == COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_Coulomb_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Coulomb_SS_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
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

                    K_Yukawa_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Yukawa_SS_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }

    /***************************************/
    /********* Regularized-Coulomb *********/
    /***************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedCoulomb_PP(num_targets_in_batch, num_sources_in_cluster,
                                batch_start, source_start,
                                target_x, target_y, target_z,
                                source_x, source_y, source_z, source_q,
                                run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedCoulomb_SS_PP(num_targets_in_batch, num_sources_in_cluster,
                                batch_start, source_start,
                                target_x, target_y, target_z, target_q,
                                source_x, source_y, source_z, source_q, source_w,
                                run_params, potential, stream_id);
                }


    /***************************************/
    /********* Regularized-Yukawa **********/
    /***************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_RegularizedYukawa_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedYukawa_SS_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                }


    /***************************************/
    /********* Atan ************************/
    /***************************************/

            } else if (run_params->kernel == ATAN) {

                K_Atan_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);



    /***************************************/
    /********* MQ **************************/
    /***************************************/

            } else if (run_params->kernel == MQ) {

                K_MQ_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);


    /***************************************/
    /********* Sin Over R ******************/
    /***************************************/

            } else if (run_params->kernel == SIN_OVER_R) {

                K_SinOverR_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

    /***************************************/
    /********* RBS U ***********************/
    /***************************************/

            } else if (run_params->kernel == RBS_U) {

                K_RBSu_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

    /***************************************/
    /********* RBS V ***********************/
    /***************************************/

            } else if (run_params->kernel == RBS_V) {

                K_RBSv_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

    /***************************************/
    /******** USER DEFINED KERNEL **********/
    /***************************************/

            } else if (run_params->kernel == USER) {

                K_User_Kernel_PP(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q,
                            run_params, potential, stream_id);

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

} /* END of Interaction_PC_Compute */
