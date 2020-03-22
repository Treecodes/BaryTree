#include <stdio.h>
#include <string.h>
#include <mpi.h>
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

#include "../comm_types/struct_comm_types.h"
#include "../comm_types/comm_types.h"

#include "../comm_windows/struct_comm_windows.h"
#include "../comm_windows/comm_windows.h"

#include "../comm_cp/comm_cp.h"

#include "../run_params/struct_run_params.h"
#include "../run_params/run_params.h"

#include "../interaction_lists/interaction_lists.h"
#include "../interaction_compute/interaction_compute.h"

#include "treedriver.h"


void treedriver(struct Particles *sources, struct Particles *targets, struct RunParams *run_params,
                double *potential, double *time_tree)
{
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    RunParams_Validate(run_params);
    Particles_ConstructOrder(sources);
    Particles_ConstructOrder(targets);
    
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    if (run_params->verbosity > 0 && rank == 0) {
        printf("[BaryTree]\n");
        printf("[BaryTree] Running BaryTree with %d ranks.\n", num_procs);
        RunParams_Print(run_params);
    }
    
    double time1;
    int total_num_direct = 0;
    int total_num_approx = 0;
    int total_num_inter = 0;
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
        STOP_TIMER(&time_tree[0]);
        
        
        START_TIMER(&time_tree[1]);
        Batches_Sources_Construct(&batches, sources, run_params);
        STOP_TIMER(&time_tree[1]);
        

        START_TIMER(&time_tree[2]);
        Clusters_Targets_Construct(&clusters, tree, run_params);
        STOP_TIMER(&time_tree[2]);
        
        
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
        STOP_TIMER(&time_tree[4]);
        
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        if (run_params->verbosity > 0) {
            total_num_approx += sum_int(local_interaction_list->num_approx, batches->numnodes);
            total_num_direct += sum_int(local_interaction_list->num_direct, batches->numnodes);
        }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

        START_TIMER(&time_tree[5]);
        InteractionCompute_CP(potential, tree, batches, local_interaction_list,
                              sources, targets, clusters, run_params);
        InteractionLists_Free(&local_interaction_list);
        STOP_TIMER(&time_tree[5]);

        //~~~~~~~~~~~~~~~~~~~~
        // Remote compute
        //~~~~~~~~~~~~~~~~~~~~
        
        if (num_procs > 1) {
                
            MPI_Barrier(MPI_COMM_WORLD);
        
            struct Tree *remote_batches = NULL;
            struct Particles *remote_sources = NULL;
            struct InteractionLists *let_interaction_list = NULL;
            
            START_TIMER(&time_tree[3]);
            Comm_CP_ConstructAndGetData(&remote_batches, &remote_sources, tree, batches, sources, run_params);
            STOP_TIMER(&time_tree[3]);
            
            START_TIMER(&time_tree[6]);
            InteractionLists_Make(&let_interaction_list, tree, remote_batches, run_params);
            STOP_TIMER(&time_tree[6]);
            
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            if (run_params->verbosity > 0) {
                total_num_approx += sum_int(let_interaction_list->num_approx, batches->numnodes);
                total_num_direct += sum_int(let_interaction_list->num_direct, batches->numnodes);
            }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            START_TIMER(&time_tree[7]);
            InteractionCompute_CP(potential, tree, remote_batches, let_interaction_list,
                                  remote_sources, targets, clusters, run_params);
            InteractionLists_Free(&let_interaction_list);
            Particles_Free(&remote_sources);
            Batches_Free(&remote_batches);
            STOP_TIMER(&time_tree[7]);
        }
        
        
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
        InteractionCompute_SubtractionPotentialCorrection(potential, targets, run_params);
        Particles_Targets_Reorder(targets, potential);
        Particles_Sources_Reorder(sources);
        STOP_TIMER(&time_tree[9]);
        
        
        //-------------------------------
        //-------------------------------
        // CLEANUP
        //-------------------------------
        //-------------------------------
        
        START_TIMER(&time_tree[10]);
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
    
        MPI_Barrier(MPI_COMM_WORLD);
        
        
        
        
//--------------------------------------------------------------------
//--------------------------------------------------------------------
// PARTICLE CLUSTER
//--------------------------------------------------------------------
//--------------------------------------------------------------------
        
    } else if (run_params->compute_type == PARTICLE_CLUSTER) {
    
        struct Tree *tree = NULL;
        struct Tree *batches = NULL;
        struct Clusters *clusters = NULL;


        //-------------------------------
        //-------------------------------
        // SETUP
        //-------------------------------
        //-------------------------------

        START_TIMER(&time_tree[0]);
        Tree_Sources_Construct(&tree, sources, run_params);
        STOP_TIMER(&time_tree[0]);
        

        START_TIMER(&time_tree[1]);
        Batches_Targets_Construct(&batches, targets, run_params);
        STOP_TIMER(&time_tree[1]);
        

        START_TIMER(&time_tree[2]);
        Clusters_Sources_Construct(&clusters, sources, tree, run_params);
        STOP_TIMER(&time_tree[2]);


        //-------------------------------
        //-------------------------------
        // COMPUTE
        //-------------------------------
        //-------------------------------
        
        struct CommTypes *comm_types = NULL;
        struct CommWindows *comm_windows = NULL;
        struct Tree **let_trees = NULL;

        struct Clusters *let_clusters = NULL;
        struct Particles *let_sources = NULL;
        
        //~~~~~~~~~~~~~~~~~~~~
        // Getting remote data
        //~~~~~~~~~~~~~~~~~~~~

        if (num_procs > 1) {
        
            MPI_Barrier(MPI_COMM_WORLD);
            
            START_TIMER(&time_tree[3]);
            CommTypesAndTrees_Construct(&comm_types, &let_trees, tree, batches, run_params);

            Particles_Alloc(&let_sources, comm_types->let_sources_length);
            Clusters_Alloc(&let_clusters, comm_types->let_clusters_length, run_params);
                                                    
            CommWindows_Create(&comm_windows, clusters, sources);
            
            for (int proc_id = 1; proc_id < num_procs; ++proc_id) {
            
                int get_from = (num_procs + rank - proc_id) % num_procs;
                
                CommWindows_Lock(comm_windows, get_from);
                //This is a non-blocking call!
                CommWindows_GetData(let_clusters, let_sources, comm_types, comm_windows, get_from, run_params);
                CommWindows_Unlock(comm_windows, get_from);
            }
            
            CommWindows_Free(&comm_windows);
            STOP_TIMER(&time_tree[3]);
        }

        //~~~~~~~~~~~~~~~~~~~~
        // Local compute
        //~~~~~~~~~~~~~~~~~~~~

        struct InteractionLists *local_interaction_list;
        
        START_TIMER(&time_tree[4]);
        InteractionLists_Make(&local_interaction_list, tree, batches, run_params);
        STOP_TIMER(&time_tree[4]);

//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        if (run_params->verbosity > 0) {
            total_num_approx += sum_int(local_interaction_list->num_approx, batches->numnodes);
            total_num_direct += sum_int(local_interaction_list->num_direct, batches->numnodes);
        }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  
        START_TIMER(&time_tree[5]);
        InteractionCompute_PC(potential, tree, batches, local_interaction_list,
                              sources, targets, clusters, run_params);
        InteractionLists_Free(&local_interaction_list);
        STOP_TIMER(&time_tree[5]);

        //~~~~~~~~~~~~~~~~~~~~
        // Remote compute
        //~~~~~~~~~~~~~~~~~~~~
        
        time_tree[6] = 0;
        time_tree[7] = 0;

        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {
        
            int get_from = (num_procs+rank-proc_id) % num_procs;
            struct InteractionLists *let_interaction_list;
            
            START_TIMER(&time1);
            InteractionLists_Make(&let_interaction_list, let_trees[get_from], batches, run_params);
            STOP_TIMER(&time1);
            time_tree[6] += time1;

//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            if (run_params->verbosity > 0) {
                total_num_approx += sum_int(let_interaction_list->num_approx, batches->numnodes);
                total_num_direct += sum_int(let_interaction_list->num_direct, batches->numnodes);
            }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            START_TIMER(&time1);
            InteractionCompute_PC(potential, let_trees[get_from], batches, let_interaction_list,
                                  let_sources, targets, let_clusters, run_params);
            InteractionLists_Free(&let_interaction_list);
            STOP_TIMER(&time1);
            time_tree[7] += time1;
            
        }
            

        //-------------------------------
        //-------------------------------
        // CORRECT AND REORDER
        //-------------------------------
        //-------------------------------
        
        time_tree[8] = 0.0;
        
        START_TIMER(&time_tree[9]);
        InteractionCompute_SubtractionPotentialCorrection(potential, targets, run_params);
        Particles_Targets_Reorder(targets, potential);
        Particles_Sources_Reorder(sources);
        STOP_TIMER(&time_tree[9]);
        

        //-------------------------------
        //-------------------------------
        // CLEANUP
        //-------------------------------
        //-------------------------------

        START_TIMER(&time_tree[10]);
        Particles_FreeOrder(sources);
        Particles_FreeOrder(targets);
        Tree_Free(&tree);
        Clusters_Free_Win(&clusters);
        Batches_Free(&batches);

        // remote pieces
        Clusters_Free(&let_clusters);

        Particles_Free(&let_sources);

        CommTypesAndTrees_Free(&comm_types, &let_trees);
        STOP_TIMER(&time_tree[10]);

        // Total setup time
        time_tree[11] = time_tree[0] + time_tree[1] + time_tree[3] + time_tree[4] + time_tree[6];
        
        // Total compute time
        time_tree[12] = time_tree[5] + time_tree[7] + time_tree[8];
    
        MPI_Barrier(MPI_COMM_WORLD);




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
        Clusters_Targets_Construct(&target_clusters, target_tree, run_params);
        STOP_TIMER(&time_tree[2]);
        
        
        //-------------------------------
        //-------------------------------
        // COMPUTE
        //-------------------------------
        //-------------------------------
    
        struct CommTypes *comm_types = NULL;
        struct CommWindows *comm_windows = NULL;
        struct Tree **let_trees = NULL;

        struct Clusters *let_clusters = NULL;
        struct Particles *let_sources = NULL;
        
        //~~~~~~~~~~~~~~~~~~~~
        // Getting remote data
        //~~~~~~~~~~~~~~~~~~~~

        if (num_procs > 1) {
        
            MPI_Barrier(MPI_COMM_WORLD);
        
            START_TIMER(&time_tree[3]);
            CommTypesAndTrees_Construct(&comm_types, &let_trees,
                                        source_tree, target_tree, run_params);

            Particles_Alloc(&let_sources, comm_types->let_sources_length);
            Clusters_Alloc(&let_clusters, comm_types->let_clusters_length, run_params);
                                                      
            CommWindows_Create(&comm_windows, source_clusters, sources);
            
            for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

                int get_from = (num_procs + rank - proc_id) % num_procs;
                
                CommWindows_Lock(comm_windows, get_from);
                //This is a non-blocking call!
                CommWindows_GetData(let_clusters, let_sources, comm_types, comm_windows, get_from, run_params);
                CommWindows_Unlock(comm_windows, get_from);
            }
            
            CommWindows_Free(&comm_windows);
            STOP_TIMER(&time_tree[3]);
        }

        //~~~~~~~~~~~~~~~~~~~~
        // Local compute
        //~~~~~~~~~~~~~~~~~~~~

        struct InteractionLists *local_interaction_list;
        
        START_TIMER(&time_tree[4]);
        InteractionLists_Make(&local_interaction_list, source_tree, target_tree, run_params);
        STOP_TIMER(&time_tree[4]);

//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        if (run_params->verbosity > 0) {
            total_num_approx += sum_int(local_interaction_list->num_approx, target_tree->numnodes);
            total_num_direct += sum_int(local_interaction_list->num_direct, target_tree->numnodes);
        }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        
        START_TIMER(&time_tree[5]);
        InteractionCompute_CC(potential, source_tree, target_tree, local_interaction_list,
                              sources, targets, source_clusters, target_clusters, run_params);
        InteractionLists_Free(&local_interaction_list);
        STOP_TIMER(&time_tree[5]);
        
        //~~~~~~~~~~~~~~~~~~~~
        // Remote compute
        //~~~~~~~~~~~~~~~~~~~~
        
        time_tree[6] = 0;
        time_tree[7] = 0;
            
        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

            int get_from = (num_procs+rank-proc_id) % num_procs;
            struct InteractionLists *let_interaction_list;
            
            START_TIMER(&time1);
            InteractionLists_Make(&let_interaction_list, let_trees[get_from], target_tree,  run_params);
            STOP_TIMER(&time1);
            time_tree[6] += time1;
             
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
            if (run_params->verbosity > 0) {
                total_num_approx += sum_int(let_interaction_list->num_approx, target_tree->numnodes);
                total_num_direct += sum_int(let_interaction_list->num_direct, target_tree->numnodes);
            }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


            START_TIMER(&time1);
            InteractionCompute_CC(potential, let_trees[get_from], target_tree, let_interaction_list,
                                  let_sources, targets, let_clusters, target_clusters, run_params);
            InteractionLists_Free(&let_interaction_list);
            STOP_TIMER(&time1);
            time_tree[7] += time1;
        }
            
            
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
        InteractionCompute_SubtractionPotentialCorrection(potential, targets, run_params);
        Particles_Targets_Reorder(targets, potential);
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
        Clusters_Free_Win(&source_clusters);
        Clusters_Free(&target_clusters);
        
        // Remote pieces
        Clusters_Free(&let_clusters);
        Particles_Free(&let_sources);
        CommTypesAndTrees_Free(&comm_types, &let_trees);
        STOP_TIMER(&time_tree[10]);
        
        // Total setup time
        time_tree[11] = time_tree[0] + time_tree[1] + time_tree[3] + time_tree[4] + time_tree[6];
        
        // Total compute time
        time_tree[12] = time_tree[5] + time_tree[7] + time_tree[8];
        
        MPI_Barrier(MPI_COMM_WORLD);
    }




//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    if (run_params->verbosity > 0) {
       
        int global_num_inter,  max_num_inter,  min_num_inter;
        int global_num_direct, max_num_direct, min_num_direct;
        int global_num_approx, max_num_approx, min_num_approx;
    
        total_num_inter = total_num_direct + total_num_approx;
        MPI_Reduce(&total_num_inter,   &global_num_inter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_inter,      &max_num_inter, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_inter,      &min_num_inter, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&total_num_direct, &global_num_direct, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_direct,    &max_num_direct, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_direct,    &min_num_direct, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&total_num_approx, &global_num_approx, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_approx,    &max_num_approx, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_approx,    &min_num_approx, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        if (rank == 0) {
            printf("[BaryTree]\n");
            printf("[BaryTree]        Cumulative interactions across all ranks: %d\n", global_num_inter);
            printf("[BaryTree]           Maximum interactions across all ranks: %d\n", max_num_inter);
            printf("[BaryTree]           Minimum interactions across all ranks: %d\n", min_num_inter);
            printf("[BaryTree]                                           Ratio: %f\n",
                   (double)max_num_inter / (double)min_num_inter);
            printf("[BaryTree]\n");
            printf("[BaryTree] Cumulative direct interactions across all ranks: %d\n", global_num_direct);
            printf("[BaryTree]    Maximum direct interactions across all ranks: %d\n", max_num_direct);
            printf("[BaryTree]    Minimum direct interactions across all ranks: %d\n", min_num_direct);
            printf("[BaryTree]                                           Ratio: %f\n",
                   (double)max_num_direct / (double)min_num_direct);
            printf("[BaryTree]\n");
            printf("[BaryTree] Cumulative approx interactions across all ranks: %d\n", global_num_approx);
            printf("[BaryTree]    Maximum approx interactions across all ranks: %d\n", max_num_approx);
            printf("[BaryTree]    Minimum approx interactions across all ranks: %d\n", min_num_approx);
            printf("[BaryTree]                                           Ratio: %f\n",
                   (double)max_num_approx / (double)min_num_approx);
            printf("[BaryTree]\n");
            printf("[BaryTree] BaryTree has finished.\n");
            printf("[BaryTree]\n");
        }
    }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
    
    return;
}
