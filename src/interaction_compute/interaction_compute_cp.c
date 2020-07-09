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

//#include "../kernels/coulomb/coulomb.h"
//#include "../kernels/yukawa/yukawa.h"
//#include "../kernels/regularized-coulomb/regularized-coulomb.h"
//#include "../kernels/regularized-yukawa/regularized-yukawa.h"
//#include "../kernels/sin-over-r/sin-over-r.h"
#include "../kernels/tcf/tcf.h"
//#include "../kernels/dcf/dcf.h"

#include "interaction_compute.h"


void InteractionCompute_CP(double *potential, struct Tree *tree, struct Tree *batches,
                           struct InteractionLists *interaction_list,
                           struct Particles *sources, struct Particles *targets,
                           struct Clusters *clusters, struct RunParams *run_params)
{
    int interp_pts_per_cluster = run_params->interp_pts_per_cluster;
    int interp_order_lim = run_params->interp_order+1;

    int num_sources   = sources->num;
    double *source_x  = sources->x;
    double *source_y  = sources->y;
    double *source_z  = sources->z;
    double *source_q  = sources->q;
    
    int num_targets   = targets->num;

    int total_num_interp_pts     = clusters->num;
    int total_num_interp_charges = clusters->num_charges;

    double *cluster_x = clusters->x;
    double *cluster_y = clusters->y;
    double *cluster_z = clusters->z;
    double *cluster_q = clusters->q;


    int **approx_inter_list = interaction_list->approx_interactions;
    int **direct_inter_list = interaction_list->direct_interactions;
    
    int *num_approx = interaction_list->num_approx;
    int *num_direct = interaction_list->num_direct;
    
    
    int *cluster_ind = tree->cluster_ind;
    
    int *target_tree_x_low_ind = tree->x_low_ind;
    int *target_tree_y_low_ind = tree->y_low_ind;
    int *target_tree_z_low_ind = tree->z_low_ind;
    
    int *target_tree_x_high_ind = tree->x_high_ind;
    int *target_tree_y_high_ind = tree->y_high_ind;
    int *target_tree_z_high_ind = tree->z_high_ind;
    
    double *target_tree_x_min = tree->x_min;
    double *target_tree_y_min = tree->y_min;
    double *target_tree_z_min = tree->z_min;
    
    
    int target_x_dim_glob = targets->xdim;
    int target_y_dim_glob = targets->ydim;
    int target_z_dim_glob = targets->zdim;
    
    double target_xdd = targets->xdd;
    double target_ydd = targets->ydd;
    double target_zdd = targets->zdd;


#ifdef OPENACC_ENABLED
    #pragma acc data copyin(source_x[0:num_sources], source_y[0:num_sources], source_z[0:num_sources], \
                            source_q[0:num_sources], \
                            cluster_x[0:total_num_interp_pts], \
                            cluster_y[0:total_num_interp_pts], \
                            cluster_z[0:total_num_interp_pts]) \
                       copy(cluster_q[0:total_num_interp_charges], \
                            potential[0:num_targets])
#endif
    {

    for (int i = 0; i < batches->numnodes; i++) {
    
        int batch_ibeg = batches->ibeg[i];
        int batch_iend = batches->iend[i];
        
        int num_approx_in_batch = num_approx[i];
        int num_direct_in_batch = num_direct[i];

        int num_sources_in_batch = batch_iend - batch_ibeg + 1;
        int batch_start =  batch_ibeg - 1;


/* * ********************************************************/
/* * ************ POTENTIAL FROM APPROX *********************/
/* * ********************************************************/

        for (int j = 0; j < num_approx_in_batch; j++) {

            int node_index = approx_inter_list[i][j];

            int cluster_charge_start = interp_pts_per_cluster*cluster_ind[node_index];
            int cluster_pts_start = interp_order_lim*cluster_ind[node_index];
            int stream_id = j%3;


    /* * *********************************************/
    /* * *************** Coulomb *********************/
    /* * *********************************************/
            if (run_params->kernel == COULOMB) {
/*
                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {
            
                        K_Coulomb_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP COULOMB SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Coulomb_CP_Hermite(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
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
*/

    /* * *********************************************/
    /* * *************** Yukawa **********************/
    /* * *********************************************/

            } else if (run_params->kernel == YUKAWA) {
/*
                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Yukawa_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);

                    } else if (run_params->singularity == SUBTRACTION) {

                        printf("**ERROR** NOT SET UP FOR CP YUKAWA SS.  EXITING.\n");
                        exit(1);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_Yukawa_CP_Hermite(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
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
*/

    /* * *************************************/
    /* * ******* Regularized Coulomb *********/
    /* * *************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {
/*
                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedCoulomb_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
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
                            source_x, source_y, source_z, source_q, source_w,
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
*/

    /* * *************************************/
    /* * ******* Regularized Yukawa **********/
    /* * *************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {
/*
                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_RegularizedYukawa_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
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
*/

    /* * *************************************/
    /* * ******* Sin Over R ******************/
    /* * *************************************/

            } else if (run_params->kernel == SIN_OVER_R) {
/*
                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_SinOverR_CP_Hermite(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);
                    }
                }
*/

    /* * *************************************/
    /* * ******* TCF *************************/
    /* * *************************************/

            } else if (run_params->kernel == TCF) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_TCF_CP_Lagrange(num_sources_in_batch, batch_start,
                            interp_pts_per_cluster, cluster_charge_start,
                            interp_order_lim, cluster_pts_start,
                            source_x, source_y, source_z, source_q,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        //K_TCF_CP_Hermite(num_sources_in_batch,
                        //    interp_pts_per_cluster, batch_start, cluster_start,
                        //    source_x, source_y, source_z, source_q,
                        //    cluster_x, cluster_y, cluster_z, cluster_q,
                        //    run_params, stream_id);
                    }
                }


    /* * *************************************/
    /* * ******* DCF *************************/
    /* * *************************************/

            } else if (run_params->kernel == DCF) {
/*
                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_DCF_CP_Lagrange(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_DCF_CP_Hermite(num_sources_in_batch,
                            interp_pts_per_cluster, batch_start, cluster_start,
                            source_x, source_y, source_z, source_q, source_w,
                            cluster_x, cluster_y, cluster_z, cluster_q,
                            run_params, stream_id);
                    }
                }
*/
            } else {
                printf("**ERROR** INVALID KERNEL. EXITING.\n");
                exit(1);
            }

        } // end loop over cluster approximations



/* * ********************************************************/
/* * ************ POTENTIAL FROM DIRECT *********************/
/* * ********************************************************/

        for (int j = 0; j < num_direct_in_batch; j++) {

            int node_index = direct_inter_list[i][j];
    
            int target_x_low_ind = target_tree_x_low_ind[node_index];
            int target_y_low_ind = target_tree_y_low_ind[node_index];
            int target_z_low_ind = target_tree_z_low_ind[node_index];
    
            int target_x_high_ind = target_tree_x_high_ind[node_index];
            int target_y_high_ind = target_tree_y_high_ind[node_index];
            int target_z_high_ind = target_tree_z_high_ind[node_index];
    
            double target_x_min = target_tree_x_min[node_index];
            double target_y_min = target_tree_y_min[node_index];
            double target_z_min = target_tree_z_min[node_index];


            int stream_id = j%3;


    /* * *********************************************/
    /* * *************** Coulomb *********************/
    /* * *********************************************/

            if (run_params->kernel == COULOMB) {
/*
                if (run_params->singularity == SKIPPING) {

                    K_Coulomb_Direct(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Coulomb_SS_Direct(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }
*/

    /* * *********************************************/
    /* * *************** Yukawa **********************/
    /* * *********************************************/

            } else if (run_params->kernel == YUKAWA) {
/*
                if (run_params->singularity == SKIPPING) {

                    K_Yukawa_Direct(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_Yukawa_SS_Direct(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }
*/

    /* * *********************************************/
    /* * ********* Regularized Coulomb ***************/
    /* * *********************************************/

            } else if (run_params->kernel == REGULARIZED_COULOMB) {
/*
                if (run_params->singularity == SKIPPING) {

                    K_RegularizedCoulomb_Direct(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedCoulomb_SS_Direct(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }
*/

    /* * *********************************************/
    /* * ********* Regularized Yukawa ****************/
    /* * *********************************************/

            } else if (run_params->kernel == REGULARIZED_YUKAWA) {
/*
                if (run_params->singularity == SKIPPING) {

                    K_RegularizedYukawa_Direct(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else if (run_params->singularity == SUBTRACTION) {

                    K_RegularizedYukawa_SS_Direct(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z, target_q,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }
*/

    /* * *********************************************/
    /* * ********* Sin Over R ************************/
    /* * *********************************************/

            } else if (run_params->kernel == SIN_OVER_R) {
/*
                if (run_params->singularity == SKIPPING) {

                    K_SinOverR_Direct(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);
                }
*/

    /* * *********************************************/
    /* * ********* TCF *******************************/
    /* * *********************************************/

            } else if (run_params->kernel == TCF) {

                if (run_params->singularity == SKIPPING) {

//                    K_TCF_Direct(num_targets_in_cluster, num_sources_in_batch,
//                            target_start, batch_start,
//                            target_x, target_y, target_z,
//                            source_x, source_y, source_z, source_q, source_w,
//                            run_params, potential, stream_id);
                            
                    K_TCF_Direct(target_x_low_ind, target_x_high_ind,
                                 target_y_low_ind, target_y_high_ind,
                                 target_z_low_ind, target_z_high_ind,
                                 target_x_min,       target_y_min,       target_z_min,

                                 target_xdd,        target_ydd,        target_zdd,
                                 target_x_dim_glob, target_y_dim_glob, target_z_dim_glob,

                                 num_sources_in_batch, batch_start,
                                 source_x, source_y, source_z, source_q,

                                 run_params, potential, stream_id);
                }


    /* * *********************************************/
    /* * ********* DCF *******************************/
    /* * *********************************************/

            } else if (run_params->kernel == DCF) {
/*
                if (run_params->singularity == SKIPPING) {

                    K_DCF_Direct(num_targets_in_cluster, num_sources_in_batch,
                            target_start, batch_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential, stream_id);
                }
*/
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
