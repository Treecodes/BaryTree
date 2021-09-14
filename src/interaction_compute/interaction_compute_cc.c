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
#include "../kernels/tcf/tcf.h"
#include "../kernels/dcf/dcf.h"

#ifdef CUDA_ENABLED
    #include "../kernels/cuda/coulomb/cuda_coulomb.h"
    #include "../kernels/cuda/tcf/cuda_tcf.h"
    #include "../kernels/cuda/dcf/cuda_dcf.h"
#endif

#include "interaction_compute.h"


void InteractionCompute_CC(double *potential, struct Tree *source_tree, struct Tree *target_tree,
                           struct InteractionLists *interaction_list,
                           struct Particles *sources, struct Particles *targets,
                           struct Clusters *source_clusters, struct Clusters *target_clusters,
                           struct RunParams *run_params)
{
    int interp_pts_per_cluster = run_params->interp_pts_per_cluster;

    int **approx_inter_list = interaction_list->approx_interactions;
    int **direct_inter_list = interaction_list->direct_interactions;
    
    int *num_approx = interaction_list->num_approx;
    int *num_direct = interaction_list->num_direct;
    
    int **source_approx_inter_list = interaction_list->cc_source_approx_interactions;
    int **target_approx_inter_list = interaction_list->cc_target_approx_interactions;
    
    int *num_source_approx = interaction_list->num_cc_source_approx;
    int *num_target_approx = interaction_list->num_cc_target_approx;
    
    
    int num_sources  = sources->num;
    double *source_x = sources->x;
    double *source_y = sources->y;
    double *source_z = sources->z;
    double *source_q = sources->q;

    int num_targets  = targets->num;


    int num_source_cluster_pts  = source_clusters->num;
    int num_source_cluster_q = source_clusters->num_charges;
    double *source_cluster_x = source_clusters->x;
    double *source_cluster_y = source_clusters->y;
    double *source_cluster_z = source_clusters->z;
    double *source_cluster_q = source_clusters->q;
    
    int num_target_cluster_pts  = target_clusters->num;
    int num_target_cluster_q = target_clusters->num_charges;
    double *target_cluster_x = target_clusters->x;
    double *target_cluster_y = target_clusters->y;
    double *target_cluster_z = target_clusters->z;
    double *target_cluster_q = target_clusters->q;
    
    
    int *source_tree_ibeg = source_tree->ibeg;
    int *source_tree_iend = source_tree->iend;
    int *source_tree_cluster_ind = source_tree->cluster_ind;
    

    int *target_tree_cluster_ind = target_tree->cluster_ind;
    
    int *target_tree_x_low_ind = target_tree->x_low_ind;
    int *target_tree_y_low_ind = target_tree->y_low_ind;
    int *target_tree_z_low_ind = target_tree->z_low_ind;
    
    int *target_tree_x_high_ind = target_tree->x_high_ind;
    int *target_tree_y_high_ind = target_tree->y_high_ind;
    int *target_tree_z_high_ind = target_tree->z_high_ind;
    
    double *target_tree_x_min = target_tree->x_min;
    double *target_tree_y_min = target_tree->y_min;
    double *target_tree_z_min = target_tree->z_min;
    
    
    int target_x_dim_glob = targets->xdim;
    int target_y_dim_glob = targets->ydim;
    int target_z_dim_glob = targets->zdim;
    
    double target_xdd = targets->xdd;
    double target_ydd = targets->ydd;
    double target_zdd = targets->zdd;
    
    
    for (int i = 0; i < target_tree->numnodes; i++) {

        int target_cluster_start = interp_pts_per_cluster * target_tree_cluster_ind[i];
        
        int num_approx_in_cluster = num_approx[i];
        int num_direct_in_cluster = num_direct[i];
        
        int num_source_approx_in_cluster = num_source_approx[i];
        int num_target_approx_in_cluster = num_target_approx[i];
        
        
        int target_x_low_ind = target_tree_x_low_ind[i];
        int target_y_low_ind = target_tree_y_low_ind[i];
        int target_z_low_ind = target_tree_z_low_ind[i];
    
        int target_x_high_ind = target_tree_x_high_ind[i];
        int target_y_high_ind = target_tree_y_high_ind[i];
        int target_z_high_ind = target_tree_z_high_ind[i];
    
        double target_x_min = target_tree_x_min[i];
        double target_y_min = target_tree_y_min[i];
        double target_z_min = target_tree_z_min[i];


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
/*
                K_Coulomb_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                    source_cluster_start, target_cluster_start,
                    source_cluster_x, source_cluster_y, source_cluster_z,
                    source_cluster_q,
                    target_cluster_x, target_cluster_y, target_cluster_z,
                    target_cluster_q,
                    run_params, stream_id);
*/


    /* * *********************************************/
    /* * ******* TCF *********************************/
    /* * *********************************************/

            } else if (run_params->kernel == TCF) {

//                K_TCF_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
//                    source_cluster_start, target_cluster_start,
//                    source_cluster_x, source_cluster_y, source_cluster_z,
//                    source_cluster_q,
//                    target_cluster_x, target_cluster_y, target_cluster_z,
//                    target_cluster_q,
//                    run_params, stream_id);
    
    
    /* * *********************************************/
    /* * ******* DCF *********************************/
    /* * *********************************************/

            } else if (run_params->kernel == DCF) {
/*
                K_DCF_CP_Lagrange(interp_pts_per_cluster, interp_pts_per_cluster,
                    source_cluster_start, target_cluster_start,
                    source_cluster_x, source_cluster_y, source_cluster_z,
                    source_cluster_q,
                    target_cluster_x, target_cluster_y, target_cluster_z,
                    target_cluster_q,
                    run_params, stream_id);
*/
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
/*
                K_Coulomb_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                    target_start, source_cluster_start,
                    target_x, target_y, target_z,
                    source_cluster_x, source_cluster_y, source_cluster_z,
                    source_cluster_q,
                    run_params, potential, stream_id);
*/


    /* * *********************************************/
    /* * ******* TCF *********************************/
    /* * *********************************************/

            } else if (run_params->kernel == TCF) {

/*
                K_TCF_PC_Lagrange(target_x_low_ind, target_x_high_ind,
                                  target_y_low_ind, target_y_high_ind,
                                  target_z_low_ind, target_z_high_ind,

                                  target_x_min,       target_y_min,     target_z_min,
                                  target_xdd,        target_ydd,        target_zdd,
                                  target_x_dim_glob, target_y_dim_glob, target_z_dim_glob,

                                  interp_pts_per_cluster, source_cluster_start,
                                  source_cluster_x, source_cluster_y, source_cluster_z, source_cluster_q,

                                  run_params, potential, stream_id);
*/
    
    
    /* * *********************************************/
    /* * ******* DCF *********************************/
    /* * *********************************************/

            } else if (run_params->kernel == DCF) {
/*
                K_DCF_PC_Lagrange(num_targets_in_cluster, interp_pts_per_cluster,
                    target_start, source_cluster_start,
                    target_x, target_y, target_z,
                    source_cluster_x, source_cluster_y, source_cluster_z,
                    source_cluster_q,
                    run_params, potential, stream_id);
*/
            } else {
                printf("**ERROR** INVALID KERNEL. EXITING.\n");
                exit(1);
            }

        } // end loop over cluster approximations



/* * ********************************************************/
/* * ************ POTENTIAL FROM TARGET APPROX (PC) *********/
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
/*
                K_Coulomb_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                    source_start, target_cluster_start,
                    source_x, source_y, source_z, source_q,
                    target_cluster_x, target_cluster_y, target_cluster_z,
                    target_cluster_q,
                    run_params, stream_id);
*/


    /* * *********************************************/
    /* * ******* TCF *********************************/
    /* * *********************************************/

            } else if (run_params->kernel == TCF) {

//                K_TCF_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
//                    source_start, target_cluster_start,
//                    source_x, source_y, source_z, source_q,
//                    target_cluster_x, target_cluster_y, target_cluster_z,
//                    target_cluster_q,
//                    run_params, stream_id);


    /* * *********************************************/
    /* * ******* DCF *********************************/
    /* * *********************************************/

            } else if (run_params->kernel == DCF) {
/*
                K_DCF_CP_Lagrange(num_sources_in_cluster, interp_pts_per_cluster,
                    source_start, target_cluster_start,
                    source_x, source_y, source_z, source_q,
                    target_cluster_x, target_cluster_y, target_cluster_z,
                    target_cluster_q,
                    run_params, stream_id);
*/
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
/*
                K_Coulomb_PP(num_targets_in_cluster, num_sources_in_cluster,
                        target_start, source_start,
                        target_x, target_y, target_z,
                        source_x, source_y, source_z, source_q,
                        run_params, potential, stream_id);
*/


    /* * *********************************************/
    /* * ********** TCF ******************************/
    /* * *********************************************/

            } else if (run_params->kernel == TCF) {

/*
                K_TCF_PP(target_x_low_ind, target_x_high_ind,
                         target_y_low_ind, target_y_high_ind,
                         target_z_low_ind, target_z_high_ind,
                         target_x_min,       target_y_min,       target_z_min,

                         target_xdd,        target_ydd,        target_zdd,
                         target_x_dim_glob, target_y_dim_glob, target_z_dim_glob,

                         num_sources_in_cluster, source_start,
                         source_x, source_y, source_z, source_q,

                         run_params, potential, stream_id);
*/
                
                
    /* * *********************************************/
    /* * ********** DCF ******************************/
    /* * *********************************************/

            } else if (run_params->kernel == DCF) {
/*
                K_DCF_PP(num_targets_in_cluster, num_sources_in_cluster,
                        target_start, source_start,
                        target_x, target_y, target_z,
                        source_x, source_y, source_z, source_q,
                        run_params, potential, stream_id);
*/
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
