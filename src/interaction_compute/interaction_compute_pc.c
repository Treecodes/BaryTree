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
#include "../kernels/tcf/tcf.h"
#include "../kernels/dcf/dcf.h"

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

    int **approx_inter_list = interaction_list->approx_interactions;
    int **direct_inter_list = interaction_list->direct_interactions;
    
    int *num_approx = interaction_list->num_approx;
    int *num_direct = interaction_list->num_direct;

    int tree_numnodes = tree->numnodes;
    int batch_numnodes = batches->numnodes;

    double *potential_direct, *potential_approx;
    make_vector(potential_direct, num_targets);
    make_vector(potential_approx, num_targets);

    memset(potential_approx, 0, num_targets * sizeof(double));
    memset(potential_direct, 0, num_targets * sizeof(double));

    int *tree_ibeg = tree->ibeg;
    int *tree_iend = tree->iend;
    int *cluster_ind = tree->cluster_ind;

#ifdef OPENACC_ENABLED
    #pragma acc data copyin(source_x[0:num_sources], source_y[0:num_sources], source_z[0:num_sources], \
                        source_q[0:num_sources], \
                        cluster_x[0:total_num_interp_pts], cluster_y[0:total_num_interp_pts], \
                        cluster_z[0:total_num_interp_pts], \
                        cluster_q[0:total_num_interp_charges]) \
                        copy(potential_approx[0:num_targets], potential_direct[0:num_targets])
#endif
    {

    for (int i = 0; i < batches->numnodes; i++) {
        int batch_ibeg = batches->ibeg[i];
        int batch_iend = batches->iend[i];

        int num_approx_in_batch = num_approx[i];
        int num_direct_in_batch = num_direct[i];

        int num_targets_in_batch = batch_iend - batch_ibeg + 1;
        int batch_start =  batch_ibeg - 1;


/* * ********************************************************/
/* * ************ POTENTIAL FROM APPROX *********************/
/* * ********************************************************/

        for (int j = 0; j < num_approx_in_batch; j++) {
            int node_index = approx_inter_list[i][j];
            int cluster_start = interp_pts_per_cluster*cluster_ind[node_index];
            int stream_id = j%3;


    /* * *********************************************/
    /* * *************** Coulomb *********************/
    /* * *********************************************/
    
            if (run_params->kernel == COULOMB) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {
            
                        K_Coulomb_PC_Lagrange(num_targets_in_batch,
                                interp_pts_per_cluster, batch_start, cluster_start,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_q,
                                run_params, potential_approx, stream_id);

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
                                run_params, potential_approx, stream_id);

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

                        K_Yukawa_PC_Lagrange(num_targets_in_batch,
                                interp_pts_per_cluster, batch_start, cluster_start,
                                target_x, target_y, target_z,
                                cluster_x, cluster_y, cluster_z, cluster_q,
                                run_params, potential_approx, stream_id);

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
                                run_params, potential_approx, stream_id);

                    } else {
                        printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                        exit(1);
                    }

                } else {
                    printf("**ERROR** INVALID CHOICE OF APPROXIMATION. EXITING. \n");
                    exit(1);
                }


    /* * *************************************/
    /* * *********** TCF *********************/
    /* * *************************************/

            } else if (run_params->kernel == TCF) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

//                        K_TCF_PC_Lagrange(num_targets_in_batch,
//                                    interp_pts_per_cluster, batch_start, cluster_start,
//                                    target_x, target_y, target_z,
//                                    cluster_x, cluster_y, cluster_z, cluster_q,
//                                    run_params, potential_approx, stream_id);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

//                        K_TCF_PC_Hermite(num_targets_in_batch,
//                                    interp_pts_per_cluster, batch_start, cluster_start,
//                                    total_num_interp_pts, target_x, target_y, target_z,
//                                    cluster_x, cluster_y, cluster_z, cluster_q,
//                                    run_params, potential_approx, stream_id);
                    }
                }


    /* * *************************************/
    /* * *********** DCF *********************/
    /* * *************************************/

            } else if (run_params->kernel == DCF) {

                if (run_params->approximation == LAGRANGE) {

                    if (run_params->singularity == SKIPPING) {

                        K_DCF_PC_Lagrange(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential_approx, stream_id);
                    }

                } else if (run_params->approximation == HERMITE) {

                    if (run_params->singularity == SKIPPING) {

                        K_DCF_PC_Hermite(num_targets_in_batch,
                                    interp_pts_per_cluster, batch_start, cluster_start,
                                    total_num_interp_pts, target_x, target_y, target_z,
                                    cluster_x, cluster_y, cluster_z, cluster_q,
                                    run_params, potential_approx, stream_id);
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

        for (int j = 0; j < num_direct_in_batch; j++) {

            int node_index = direct_inter_list[i][j];
            int source_start = tree_ibeg[node_index] - 1;
            int source_end = tree_iend[node_index];
            int num_sources_in_cluster = source_end-source_start;
            int stream_id = j%3;


    /* * *********************************************/
    /* * *************** Coulomb *********************/
    /* * *********************************************/

            if (run_params->kernel == COULOMB) {

                if (run_params->singularity == SKIPPING) {

                    K_Coulomb_Direct(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential_direct, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }


    /* * *********************************************/
    /* * *************** Yukawa **********************/
    /* * *********************************************/

            } else if (run_params->kernel == YUKAWA) {

                if (run_params->singularity == SKIPPING) {

                    K_Yukawa_Direct(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential_direct, stream_id);

                } else {
                    printf("**ERROR** INVALID CHOICE OF SINGULARITY. EXITING. \n");
                    exit(1);
                }


    /* * *************************************/
    /* * ******* TCF *************************/
    /* * *************************************/

            } else if (run_params->kernel == TCF) {

//                K_TCF_Direct(num_targets_in_batch, num_sources_in_cluster,
//                            batch_start, source_start,
//                            target_x, target_y, target_z,
//                            source_x, source_y, source_z, source_q, source_w,
//                            run_params, potential_direct, stream_id);


    /* * *************************************/
    /* * ******* DCF *************************/
    /* * *************************************/

            } else if (run_params->kernel == DCF) {

                K_DCF_Direct(num_targets_in_batch, num_sources_in_cluster,
                            batch_start, source_start,
                            target_x, target_y, target_z,
                            source_x, source_y, source_z, source_q, source_w,
                            run_params, potential_direct, stream_id);

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

    for (int k = 0; k < num_targets; k++) 
        potential[k] += potential_direct[k];

    for (int k = 0; k < num_targets; k++) 
        potential[k] += potential_approx[k];

    free_vector(potential_direct);
    free_vector(potential_approx);

    return;

} /* END of Interaction_PC_Compute */
