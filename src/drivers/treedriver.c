#include <stdio.h>
#include <string.h>
#include <limits.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../utilities/timers.h"
#include "../utilities/enums.h"

#include "../tree/struct_tree.h"
#include "../tree/tree.h"
#include "../tree/batches.h"

#include "../particles/struct_particles.h"
#include "../particles/particles.h"

#include "../clusters/struct_clusters.h"
#include "../clusters/clusters.h"

#include "../run_params/struct_run_params.h"
#include "../run_params/run_params.h"

#include "../interaction_lists/interaction_lists.h"
#include "../interaction_compute/interaction_compute.h"

#include "treedriver.h"


void treedriver(struct Particles *sources, struct Particles *targets, struct RunParams *run_params,
                double *potential, double *time_tree)
{
    RunParams_Validate(run_params);
    Particles_ConstructOrder(sources);
    Particles_ConstructOrder(targets);
    
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    if (run_params->verbosity > 0) {
        printf("[BaryTree]\n");
        printf("[BaryTree] Running BaryTree.\n");
        RunParams_Print(run_params);
    }
    
    long long int total_num_direct = 0;
    long long int total_num_approx = 0;
    long long int total_num_inter = 0;

    long long int total_num_direct_interact = 0;
    long long int total_num_approx_interact = 0;
    long long int total_num_interact = 0;

    // These types of interactions only occur for CC
    long long int total_num_source_approx = 0;
    long long int total_num_target_approx = 0;
    long long int total_num_source_approx_interact = 0;
    long long int total_num_target_approx_interact = 0;
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~




//--------------------------------------------------------------------
//--------------------------------------------------------------------
// CLUSTER PARTICLE
//--------------------------------------------------------------------
//--------------------------------------------------------------------
    
    if (run_params->compute_type == CLUSTER_PARTICLE) {
    
        struct Tree *tree = NULL;
        struct Tree *batches = NULL;
        struct Clusters *clusters = NULL;

        
        //-------------------------------
        //-------------------------------
        // SETUP
        //-------------------------------
        //-------------------------------


        START_TIMER(&time_tree[0]);
        Tree_Targets_Construct(&tree, targets, run_params);
#ifdef OPENACC_ENABLED
        #pragma acc enter data create(potential[0:targets->num])
#endif
        STOP_TIMER(&time_tree[0]);
        
        START_TIMER(&time_tree[1]);
        Batches_Sources_Construct(&batches, sources, run_params);
#ifdef OPENACC_ENABLED
        #pragma acc enter data copyin(sources->x[0:sources->num], sources->y[0:sources->num], \
                                      sources->z[0:sources->num], sources->q[0:sources->num])
#endif
        STOP_TIMER(&time_tree[1]);

//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        if (run_params->verbosity > 0) {
            Tree_Print(tree);
            Batches_Print(batches);
        }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        
        
        //-------------------------------
        //-------------------------------
        // COMPUTE
        //-------------------------------
        //-------------------------------
        
        //~~~~~~~~~~~~~~~~~~~~
        // Local compute
        //~~~~~~~~~~~~~~~~~~~~
        
        struct InteractionLists *local_interaction_list = NULL;

        START_TIMER(&time_tree[4]);
        InteractionLists_Make(&local_interaction_list, tree, batches, run_params);
        Tree_Downpass_Interact(tree);
        STOP_TIMER(&time_tree[4]);

        START_TIMER(&time_tree[2]);
        Clusters_Targets_Construct(&clusters, tree, run_params);
        STOP_TIMER(&time_tree[2]);
        
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        if (run_params->verbosity > 0) {
            total_num_approx += sum_int(local_interaction_list->num_approx, batches->numnodes);
            total_num_direct += sum_int(local_interaction_list->num_direct, batches->numnodes);

            for (int i = 0; i < batches->numnodes; ++i) {
                for (int j = 0; j < local_interaction_list->num_direct[i]; ++j) {
                    total_num_direct_interact += (long long int) batches->numpar[i]
                        * (long long int) tree->numpar[local_interaction_list->direct_interactions[i][j]];
                }
                for (int j = 0; j < local_interaction_list->num_approx[i]; ++j) {
                    total_num_approx_interact += (long long int) batches->numpar[i]
                        * (long long int) run_params->interp_pts_per_cluster;
                }
            }
        }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        START_TIMER(&time_tree[5]);
        InteractionCompute_CP(potential, tree, batches, local_interaction_list,
                              sources, targets, clusters, run_params);
        InteractionLists_Free(&local_interaction_list);
        STOP_TIMER(&time_tree[5]);


        //-------------------------------
        //-------------------------------
        // DOWNPASS
        //-------------------------------
        //-------------------------------

        START_TIMER(&time_tree[8]);
        InteractionCompute_Downpass(potential, tree, targets, clusters, run_params);
        STOP_TIMER(&time_tree[8]);


        //-------------------------------
        //-------------------------------
        // CORRECT AND REORDER
        //-------------------------------
        //-------------------------------

        START_TIMER(&time_tree[9]);
        Particles_Sources_Reorder(sources);
        STOP_TIMER(&time_tree[9]);

        
        //-------------------------------
        //-------------------------------
        // CLEANUP
        //-------------------------------
        //-------------------------------
        
        START_TIMER(&time_tree[10]);
#ifdef OPENACC_ENABLED
        #pragma acc exit data copyout(potential[0:targets->num])
        #pragma acc exit data delete(sources->x, sources->y, sources->z, sources->q)
#endif
        Particles_FreeOrder(sources);
        Particles_FreeOrder(targets);
        Tree_Free(&tree);
        Clusters_Free(&clusters);
        Batches_Free(&batches);
        STOP_TIMER(&time_tree[10]);
        
        // Total setup time
        time_tree[11] = time_tree[0] + time_tree[1] + time_tree[3] + time_tree[4] + time_tree[6];
        
        // Total compute time
        time_tree[12] = time_tree[5] + time_tree[7] + time_tree[8];
    
        
        
        
//--------------------------------------------------------------------
//--------------------------------------------------------------------
//CLUSTER CLUSTER
//--------------------------------------------------------------------
//--------------------------------------------------------------------

    } else if (run_params->compute_type == CLUSTER_CLUSTER) {
    
        struct Tree *source_tree = NULL;
        struct Tree *target_tree = NULL;

        struct Clusters *source_clusters = NULL;
        struct Clusters *target_clusters = NULL;


        //-------------------------------
        //-------------------------------
        // SETUP
        //-------------------------------
        //-------------------------------

        START_TIMER(&time_tree[0]);
        Tree_Sources_Construct(&source_tree, sources, run_params);
        STOP_TIMER(&time_tree[0]);

        START_TIMER(&time_tree[1]);
        Tree_Targets_Construct(&target_tree, targets, run_params);
        STOP_TIMER(&time_tree[1]);
         
        START_TIMER(&time_tree[2]);
        Clusters_Sources_Construct(&source_clusters, sources, source_tree, run_params);
        STOP_TIMER(&time_tree[2]);
        
        //~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        if (run_params->verbosity > 0) {
            Tree_Print(source_tree);
            Tree_Print(target_tree);
        }
        //~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        
        
        //-------------------------------
        //-------------------------------
        // COMPUTE
        //-------------------------------
        //-------------------------------
    
        //~~~~~~~~~~~~~~~~~~~~
        // Local compute
        //~~~~~~~~~~~~~~~~~~~~

        struct InteractionLists *local_interaction_list;
        
        START_TIMER(&time_tree[4]);
        InteractionLists_Make(&local_interaction_list, source_tree, target_tree, run_params);
        STOP_TIMER(&time_tree[4]);

        double temp_time;
        START_TIMER(&temp_time);
        Clusters_Targets_Construct(&target_clusters, target_tree, run_params);
        STOP_TIMER(&temp_time);
        time_tree[2] += temp_time;

//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        if (run_params->verbosity > 0) {
            total_num_approx += sum_int(local_interaction_list->num_approx, target_tree->numnodes);
            total_num_direct += sum_int(local_interaction_list->num_direct, target_tree->numnodes);
            
            total_num_source_approx += sum_int(local_interaction_list->num_cc_source_approx,
                                               target_tree->numnodes);
            total_num_target_approx += sum_int(local_interaction_list->num_cc_target_approx,
                                               target_tree->numnodes);


            total_num_approx_interact += (long long int) sum_int(local_interaction_list->num_approx, target_tree->numnodes)
                                       * (long long int) run_params->interp_pts_per_cluster * run_params->interp_pts_per_cluster;
  
            for (int i = 0; i < target_tree->numnodes; ++i) {
                for (int j = 0; j < local_interaction_list->num_direct[i]; ++j) {
                    total_num_direct_interact += (long long int) target_tree->numpar[i]
                        * (long long int) source_tree->numpar[local_interaction_list->direct_interactions[i][j]];
                }
                for (int j = 0; j < local_interaction_list->num_cc_source_approx[i]; ++j) {
                    total_num_source_approx_interact += (long long int) target_tree->numpar[i]
                        * (long long int) run_params->interp_pts_per_cluster;
                }
                for (int j = 0; j < local_interaction_list->num_cc_target_approx[i]; ++j) {
                    total_num_target_approx_interact += (long long int) run_params->interp_pts_per_cluster
                        * (long long int) source_tree->numpar[local_interaction_list->cc_target_approx_interactions[i][j]];
                }
            }
        }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        
        START_TIMER(&time_tree[5]);
        InteractionCompute_CC(potential, source_tree, target_tree, local_interaction_list,
                              sources, targets, source_clusters, target_clusters, run_params);
        InteractionLists_Free(&local_interaction_list);
        STOP_TIMER(&time_tree[5]);
        
            
        //-------------------------------
        //-------------------------------
        // DOWNPASS
        //-------------------------------
        //-------------------------------
        
        START_TIMER(&time_tree[8]);
        InteractionCompute_Downpass(potential, target_tree, targets, target_clusters, run_params);
        STOP_TIMER(&time_tree[8]);
            
        
        //-------------------------------
        //-------------------------------
        //CORRECT AND REORDER
        //-------------------------------
        //-------------------------------

        START_TIMER(&time_tree[9]);
        Particles_Sources_Reorder(sources);
        STOP_TIMER(&time_tree[9]);
               
               
        //-------------------------------
        //-------------------------------
        //CLEANUP
        //-------------------------------
        //-------------------------------

        START_TIMER(&time_tree[10]);
        Particles_FreeOrder(sources);
        Particles_FreeOrder(targets);
        Tree_Free(&source_tree);
        Tree_Free(&target_tree);
        Clusters_Free(&source_clusters);
        Clusters_Free(&target_clusters);
        STOP_TIMER(&time_tree[10]);
        
        // Total setup time
        time_tree[11] = time_tree[0] + time_tree[1] + time_tree[3] + time_tree[4] + time_tree[6];
        
        // Total compute time
        time_tree[12] = time_tree[5] + time_tree[7] + time_tree[8];
    }




//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    if (run_params->verbosity > 0) {
       
        total_num_inter = total_num_direct + total_num_approx
                        + total_num_source_approx + total_num_target_approx;
                        
        printf("[BaryTree]\n");
        printf("[BaryTree] Interaction information: \n");
        printf("[BaryTree]\n");
        printf("[BaryTree]        Cumulative interactions: %d\n", total_num_inter);
        printf("[BaryTree]\n");
        printf("[BaryTree] Cumulative direct interactions: %d\n", total_num_direct);
        printf("[BaryTree]\n");
        printf("[BaryTree] Cumulative approx interactions: %d\n", total_num_approx);
        printf("[BaryTree]\n");
        
        // These types of interactions only occur for CC
        if (run_params->compute_type == CLUSTER_CLUSTER) {
            printf("[BaryTree] Cumulative source approx inter: %d\n", total_num_source_approx);
            printf("[BaryTree]\n");
            printf("[BaryTree] Cumulative target approx inter: %d\n", total_num_target_approx);
            printf("[BaryTree]\n");
        }


        /* For the pointwise interactions */

        total_num_interact = total_num_direct_interact + total_num_approx_interact
                           + total_num_source_approx_interact + total_num_target_approx_interact;

        printf("[BaryTree]\n");
        printf("[BaryTree]               Cumulative pointwise interactions: %lld\n", total_num_interact);
        printf("[BaryTree]\n"); 

        printf("[BaryTree]        Cumulative direct pointwise interactions: %lld\n", total_num_direct_interact);
        printf("[BaryTree]\n");

        printf("[BaryTree]        Cumulative approx pointwise interactions: %lld\n", total_num_approx_interact);
        printf("[BaryTree]\n"); 

        // These types of interactions only occur for CC
        if (run_params->compute_type == CLUSTER_CLUSTER) {
            printf("[BaryTree] Cumulative source approx pointwise interactions: %lld\n", total_num_source_approx_interact);
            printf("[BaryTree]\n");

            printf("[BaryTree] Cumulative target approx pointwise interactions: %lld\n", total_num_target_approx_interact);
            printf("[BaryTree]\n");
        }
        
        printf("[BaryTree] BaryTree has finished.\n");
        printf("[BaryTree]\n");
    }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
    return;
}
