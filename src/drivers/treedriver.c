#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <limits.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../utilities/timers.h"
#include "../utilities/enums.h"

#include "../globvars.h"

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

//Just for the moment to put math stuff
#include <math.h>

double *tt, *ww;


void treedriver(struct Particles *sources, struct Particles *targets, struct RunParams *run_params,
                double *potential, double *time_tree)
{


    /* changing values of our extern variables */
    make_vector(tt, run_params->interp_order+1);
    make_vector(ww, run_params->interp_order+1);

    for (int i = 0; i < run_params->interp_order+1; i++)
        tt[i] = cos(i * M_PI / run_params->interp_order);

    ww[0] = 0.25 * (run_params->interp_order*run_params->interp_order/3.0 + 1.0/6.0);
    ww[run_params->interp_order] = -ww[0];

    for (int i = 1; i < run_params->interp_order; i++) {
        double xx = i * M_PI / run_params->interp_order;
        ww[i] = -cos(xx) / (2 * sin(xx) * sin(xx));
    }
    
    
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    RunParams_Validate(run_params);
    Particles_ConstructOrder(sources);
    Particles_ConstructOrder(targets);
    
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
    if (run_params->verbosity > 0 && rank == 0) {
        printf("\nRunning BaryTree with %d ranks.\n", num_procs);
        RunParams_Print(run_params);
    }
    
    double time1;
    int total_num_direct = 0;
    int total_num_approx = 0;
    int total_num_interactions = 0;
    int global_num_interactions = 0;
    int max_num_interactions = 0;
    int min_num_interactions = 0;
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@

    

    //--------------------------------------
    //--------------------------------------
    //CLUSTER PARTICLE
    //--------------------------------------
    //--------------------------------------
    

    /* call setup to allocate arrays for Taylor expansions and setup global vars */
    if (run_params->compute_type == CLUSTER_PARTICLE) {
    
        struct Tree *tree = NULL;
        struct Tree *batches = NULL;
        struct Clusters *clusters = NULL;
        

        START_TIMER(&time_tree[0]);
        Tree_Targets_Construct(&tree, targets, run_params);
        STOP_TIMER(&time_tree[0]);
        
        
        START_TIMER(&time_tree[1]);
        Batches_Sources_Construct(&batches, sources, run_params);
        STOP_TIMER(&time_tree[1]);
        

        START_TIMER(&time_tree[2]);
        Clusters_Targets_Construct(&clusters, tree, run_params);
        STOP_TIMER(&time_tree[2]);
        
        
        //-------------------
        //BEGIN COMPUTE PHASE
        //-------------------
        
        struct InteractionLists *local_interaction_list = NULL;

        // Beginning local computation

        START_TIMER(&time_tree[4]);
        InteractionLists_Make(&local_interaction_list, tree, batches, run_params);
        STOP_TIMER(&time_tree[4]);
        
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
        if (run_params->verbosity > 0) {
            total_num_approx += sum_int(local_interaction_list->num_approx, batches->numnodes);
            total_num_direct += sum_int(local_interaction_list->num_direct, batches->numnodes);
        }
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@

        START_TIMER(&time_tree[5]);
        InteractionCompute_CP(tree, batches,
                        local_interaction_list,
                        sources->x, sources->y, sources->z, sources->q, sources->w,
                        targets->x, targets->y, targets->z, targets->q,
                        clusters->x, clusters->y, clusters->z, clusters->q, clusters->w,
                        potential, sources->num, targets->num, clusters->num,
                        run_params);
        InteractionLists_Free(local_interaction_list);
        STOP_TIMER(&time_tree[5]);


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
            
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
            if (run_params->verbosity > 0) {
                total_num_approx += sum_int(local_interaction_list->num_approx, batches->numnodes);
                total_num_direct += sum_int(local_interaction_list->num_direct, batches->numnodes);
            }
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@

            START_TIMER(&time_tree[7]);
            InteractionCompute_CP(tree, remote_batches,
                                  let_interaction_list,
                                  remote_sources->x, remote_sources->y, remote_sources->z, remote_sources->q, remote_sources->w,
                                  targets->x, targets->y, targets->z, targets->q,
                                  clusters->x, clusters->y, clusters->z, clusters->q, clusters->w,
                                  potential, remote_sources->num, targets->num, clusters->num,
                                  run_params);
                                  
            InteractionLists_Free(let_interaction_list);
            Particles_Free(remote_sources);
            Batches_Free(remote_batches);
            STOP_TIMER(&time_tree[7]);
        }


        START_TIMER(&time_tree[8]);
        InteractionCompute_Downpass(tree,
                                targets->x, targets->y, targets->z, targets->q,
                                clusters->x, clusters->y, clusters->z, clusters->q, clusters->w,
                                potential,
                                targets->num, clusters->num_charges, clusters->num_weights,
                                run_params);
        STOP_TIMER(&time_tree[8]);


        START_TIMER(&time_tree[9]);
        InteractionCompute_SubtractionPotentialCorrection(potential, targets->q, targets->num, run_params);
        Particles_Targets_Reorder(targets, potential);
        Particles_Sources_Reorder(sources);
        STOP_TIMER(&time_tree[9]);
        
        
        //-------------------
        //CLEANUP
        //-------------------
        
        START_TIMER(&time_tree[10]);
        Particles_FreeOrder(sources);
        Particles_FreeOrder(targets);
        Tree_Free(tree);
        Clusters_Free(clusters);
        Batches_Free(batches);
        STOP_TIMER(&time_tree[10]);
        
        // Total setup time
        time_tree[11] = time_tree[0] + time_tree[1] + time_tree[3] + time_tree[4] + time_tree[6];
        
        // Total compute time
        time_tree[12] = time_tree[5] + time_tree[7] + time_tree[8];
    
        MPI_Barrier(MPI_COMM_WORLD);

        return;
        
        
        
        //--------------------------------------
        //--------------------------------------
        //PARTICLE CLUSTER
        //--------------------------------------
        //--------------------------------------
        
    } else if (run_params->compute_type == PARTICLE_CLUSTER) {
    
        struct Tree *tree = NULL;
        struct Tree *batches = NULL;
        struct Clusters *clusters = NULL;


        START_TIMER(&time_tree[0]);
        Tree_Sources_Construct(&tree, sources, run_params);
        STOP_TIMER(&time_tree[0]);
        

        START_TIMER(&time_tree[1]);
        Batches_Targets_Construct(&batches, targets, run_params);
        STOP_TIMER(&time_tree[1]);
        

        START_TIMER(&time_tree[2]);
        Clusters_Sources_Construct(&clusters, sources, tree, run_params);
        STOP_TIMER(&time_tree[2]);


        //-------------------
        //BEGIN COMPUTE PHASE
        //-------------------

        MPI_Barrier(MPI_COMM_WORLD);
        
        struct CommTypes *comm_types = NULL;
        struct CommWindows *comm_windows = NULL;
        struct Tree **let_trees = NULL;

        struct Clusters *let_clusters = NULL;
        struct Particles *let_sources = NULL;
        
        
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
        
        CommWindows_Free(comm_windows);
        STOP_TIMER(&time_tree[3]);


        //-------------------
        //LOCAL COMPUTE
        //-------------------
        
        struct InteractionLists *local_interaction_list;
        
        START_TIMER(&time_tree[4]);
        InteractionLists_Make(&local_interaction_list, tree, batches, run_params);
        STOP_TIMER(&time_tree[4]);

//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
        if (run_params->verbosity > 0) {
            total_num_approx += sum_int(local_interaction_list->num_approx, batches->numnodes);
            total_num_direct += sum_int(local_interaction_list->num_direct, batches->numnodes);
        }
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
        
        START_TIMER(&time_tree[4]);
        InteractionCompute_PC(tree, batches,
                        local_interaction_list,
                        sources->x, sources->y, sources->z, sources->q, sources->w,
                        targets->x, targets->y, targets->z, targets->q,
                        clusters->x, clusters->y, clusters->z, clusters->q, clusters->w,
                        potential,
                        sources->num, targets->num, clusters->num,
                        run_params);
        InteractionLists_Free(local_interaction_list);
        STOP_TIMER(&time_tree[4]);


        //-------------------
        //REMOTE COMPUTE
        //-------------------
        
        time_tree[6] = 0;
        time_tree[7] = 0;
            
        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {
        
            int get_from = (num_procs+rank-proc_id) % num_procs;
            struct InteractionLists *let_interaction_list;
            
            START_TIMER(&time1);
            InteractionLists_Make(&let_interaction_list, let_trees[get_from], batches, run_params);
            STOP_TIMER(&time1);
            time_tree[6] += time1;

//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
            if (run_params->verbosity > 0) {
                total_num_approx += sum_int(local_interaction_list->num_approx, batches->numnodes);
                total_num_direct += sum_int(local_interaction_list->num_direct, batches->numnodes);
            }
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@

            START_TIMER(&time1);
            InteractionCompute_PC(let_trees[get_from], batches,
                                   let_interaction_list,
                                   let_sources->x, let_sources->y, let_sources->z, let_sources->q, let_sources->w,
                                   targets->x, targets->y, targets->z, targets->q,
                                   let_clusters->x, let_clusters->y, let_clusters->z, let_clusters->q, let_clusters->w,
                                   potential,
                                   let_sources->num, targets->num, let_clusters->num,
                                   run_params);
            InteractionLists_Free(let_interaction_list);
            STOP_TIMER(&time1);
            time_tree[7] += time1;
            
        }
            
        
        //-------------------
        //CORRECT AND REORDER
        //-------------------
        
        time_tree[8] = 0.0;
        
        START_TIMER(&time_tree[9]);
        InteractionCompute_SubtractionPotentialCorrection(potential, targets->q, targets->num, run_params);
        Particles_Targets_Reorder(targets, potential);
        Particles_Sources_Reorder(sources);
        STOP_TIMER(&time_tree[9]);
        
        
        //-------------------
        //CLEANUP
        //-------------------

        START_TIMER(&time_tree[10]);
        Particles_FreeOrder(sources);
        Particles_FreeOrder(targets);
        Tree_Free(tree);
        Clusters_Free_Win(clusters);
        Batches_Free(batches);
        
        //remote pieces
        Clusters_Free(let_clusters);
        Particles_Free(let_sources);
        CommTypesAndTrees_Free(comm_types, let_trees);
        STOP_TIMER(&time_tree[10]);
        
        // Total setup time
        time_tree[11] = time_tree[0] + time_tree[1] + time_tree[3] + time_tree[4] + time_tree[6];
        
        //Total compute time
        time_tree[12] = time_tree[5] + time_tree[7] + time_tree[8];
    
        MPI_Barrier(MPI_COMM_WORLD);

        return;



        //--------------------------------------
        //--------------------------------------
        //CLUSTER CLUSTER
        //--------------------------------------
        //--------------------------------------

    } else if (run_params->compute_type == CLUSTER_CLUSTER) {
    
        struct Tree *source_tree = NULL;
        struct Tree *target_tree = NULL;

        struct Clusters *source_clusters = NULL;
        struct Clusters *target_clusters = NULL;


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
        
        
        //-------------------
        //COMPUTE PHASE
        //-------------------
        
        MPI_Barrier(MPI_COMM_WORLD);
    
        struct CommTypes *comm_types = NULL;
        struct CommWindows *comm_windows = NULL;
        struct Tree **let_trees = NULL;

        struct Clusters *let_clusters = NULL;
        struct Particles *let_sources = NULL;


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
        
        CommWindows_Free(comm_windows);
        STOP_TIMER(&time_tree[3]);

        
        //-------------------
        //LOCAL COMPUTE
        //-------------------

        struct InteractionLists *local_interaction_list;
        
        START_TIMER(&time_tree[4]);
        InteractionLists_Make(&local_interaction_list, source_tree, target_tree, run_params);
        STOP_TIMER(&time_tree[4]);

//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
        if (run_params->verbosity > 0) {
            total_num_approx += sum_int(local_interaction_list->num_approx, target_tree->numnodes);
            total_num_direct += sum_int(local_interaction_list->num_direct, target_tree->numnodes);
        }
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
        
        START_TIMER(&time_tree[5]);
        InteractionCompute_CC(source_tree, target_tree,
                        local_interaction_list,
                        sources->x, sources->y, sources->z, sources->q, sources->w,
                        targets->x, targets->y, targets->z, targets->q,
                        source_clusters->x, source_clusters->y, source_clusters->z,
                        source_clusters->q, source_clusters->w,
                        target_clusters->x, target_clusters->y, target_clusters->z,
                        target_clusters->q, target_clusters->w,
                        potential,
                        sources->num, targets->num, source_clusters->num, target_clusters->num,
                        run_params);
        InteractionLists_Free(local_interaction_list);
        STOP_TIMER(&time_tree[5]);
        

        //-------------------
        //REMOTE COMPUTE
        //-------------------
        
        time_tree[6] = 0;
        time_tree[7] = 0;
            
        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

            int get_from = (num_procs+rank-proc_id) % num_procs;
            struct InteractionLists *let_interaction_list;
            
            START_TIMER(&time1);
            InteractionLists_Make(&let_interaction_list, let_trees[get_from], target_tree,  run_params);
            STOP_TIMER(&time1);
            time_tree[6] += time1;
             
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@
            if (run_params->verbosity > 0) {
                total_num_approx += sum_int(local_interaction_list->num_approx, target_tree->numnodes);
                total_num_direct += sum_int(local_interaction_list->num_direct, target_tree->numnodes);
            }
//~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@~@

            START_TIMER(&time1);
            InteractionCompute_CC(let_trees[get_from], target_tree,
                                   let_interaction_list,
                                   let_sources->x, let_sources->y, let_sources->z,
                                   let_sources->q, let_sources->w,
                                   targets->x, targets->y, targets->z, targets->q,
                                   let_clusters->x, let_clusters->y, let_clusters->z,
                                   let_clusters->q, let_clusters->w,
                                   target_clusters->x, target_clusters->y, target_clusters->z,
                                   target_clusters->q, target_clusters->w,
                                   potential,
                                   let_sources->num, targets->num, let_clusters->num, target_clusters->num,
                                   run_params);
            InteractionLists_Free(let_interaction_list);
            STOP_TIMER(&time1);
            time_tree[7] += time1;
        }
            
        
        //-------------------
        //DOWNPASS
        //-------------------
        
        START_TIMER(&time_tree[8]);
        InteractionCompute_Downpass(target_tree,
                                targets->x, targets->y, targets->z, targets->q,
                                target_clusters->x, target_clusters->y, target_clusters->z,
                                target_clusters->q, target_clusters->w,
                                potential,
                                targets->num, target_clusters->num_charges, target_clusters->num_weights,
                                run_params);
        STOP_TIMER(&time_tree[8]);
            
        
        //-------------------
        //CORRECT AND REORDER
        //-------------------

        START_TIMER(&time_tree[9]);
        InteractionCompute_SubtractionPotentialCorrection(potential, targets->q, targets->num, run_params);
        Particles_Targets_Reorder(targets, potential);
        Particles_Sources_Reorder(sources);
        STOP_TIMER(&time_tree[9]);
               
               
        //-------------------
        //CLEANUP
        //-------------------

        START_TIMER(&time_tree[10]);
        Particles_FreeOrder(sources);
        Particles_FreeOrder(targets);
        Tree_Free(source_tree);
        Tree_Free(target_tree);
        Clusters_Free_Win(source_clusters);
        Clusters_Free(target_clusters);
        
        //remote pieces
        Clusters_Free(let_clusters);
        Particles_Free(let_sources);
        CommTypesAndTrees_Free(comm_types, let_trees);
        STOP_TIMER(&time_tree[10]);
        
        //Total setup time
        time_tree[11] = time_tree[0] + time_tree[1] + time_tree[3] + time_tree[4] + time_tree[6];
        
        //Total compute time
        time_tree[12] = time_tree[5] + time_tree[7] + time_tree[8];
        
        MPI_Barrier(MPI_COMM_WORLD);

        return;
    }
    
    
    

//    if (run_params->verbosity > 0) {
//        printf("Tree information: \n\n");
//
//        printf("                      numpar: %d\n", troot->numpar);
//        printf("                       x_mid: %e\n", troot->x_mid);
//        printf("                       y_mid: %e\n", troot->y_mid);
//        printf("                       z_mid: %e\n\n", troot->z_mid);
//        printf("                      radius: %f\n\n", troot->radius);
//        printf("                       x_len: %e\n", troot->x_max - troot->x_min);
//        printf("                       y_len: %e\n", troot->y_max - troot->y_min);
//        printf("                       z_len: %e\n\n", troot->z_max - troot->z_min);
//        printf("                      torder: %d\n", interpolationOrder);
//        printf("                       theta: %f\n", theta);
//        printf("                  maxparnode: %d\n", maxparnode);
//        printf("            number of leaves: %d\n", numleaves);
//        printf("             number of nodes: %d\n", numnodes);
//        printf("           target batch size: %d\n", batch_size);
//        printf("           number of batches: %d\n\n", batches->numnodes);
//    }




//    if (run_params->verbosity > 0) {
//        total_num_interactions=total_num_direct+total_num_approx;
//        printf("Interaction information: \n");
//        printf("rank %d: number of direct batch-cluster interactions: %d\n", rank, total_num_approx);
//        printf("rank %d: number of approx batch-cluster interactions: %d\n", rank, total_num_direct);
//        printf("rank %d:  total number of batch-cluster interactions: %d\n\n", rank, total_num_interactions);
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//
//    if (run_params -> verbosity > 0) {
//        MPI_Reduce(&total_num_interactions,&global_num_interactions, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//        MPI_Reduce(&total_num_interactions,&max_num_interactions, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
//        MPI_Reduce(&total_num_interactions,&min_num_interactions, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
//        if (rank==0){
//            printf("Cumulative number of interactions across all ranks: %d\n", global_num_interactions);
//            printf("   Maximum number of interactions across all ranks: %d\n", max_num_interactions);
//            printf("   Minimum number of interactions across all ranks: %d\n", min_num_interactions);
//            printf("                                             Ratio: %f\n\n", (double)max_num_interactions/(double)min_num_interactions );
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }

} /* END function treecode */
