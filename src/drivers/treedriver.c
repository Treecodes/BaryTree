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
    Particles_Validate(sources, targets);
    Particles_ConstructOrder(sources);
    Particles_ConstructOrder(targets);

    int sources_w_dummy = 0; 
    int targets_q_dummy = 0; 

    if (sources->w == NULL) {
        sources_w_dummy = 1; 
        make_vector(sources->w, sources->num);
    }    
    if (targets->q == NULL) {
        targets_q_dummy = 1; 
        make_vector(targets->q, targets->num);
    }  
    
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    if (run_params->verbosity > 0 && rank == 0) {
        printf("[BaryTree]\n");
        printf("[BaryTree] Running BaryTree with %d ranks.\n", num_procs);
        RunParams_Print(run_params);
    }
    
    double time1;
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
        #pragma acc enter data copyin(targets->x[0:targets->num], targets->y[0:targets->num], \
                                      targets->z[0:targets->num])
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc enter data copyin(targets->q[0:targets->num])
        }
        #pragma acc enter data create(potential[0:targets->num])
#endif
        STOP_TIMER(&time_tree[0]);
        
        START_TIMER(&time_tree[1]);
        Batches_Sources_Construct(&batches, sources, run_params);
#ifdef OPENACC_ENABLED
        #pragma acc enter data copyin(sources->x[0:sources->num], sources->y[0:sources->num], \
                                      sources->z[0:sources->num], sources->q[0:sources->num])
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc enter data copyin(sources->w[0:sources->num])
        }
#endif
        STOP_TIMER(&time_tree[1]);

        START_TIMER(&time_tree[2]);
        Clusters_Targets_Construct(&clusters, targets, tree, run_params);
        STOP_TIMER(&time_tree[2]);
        
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
        STOP_TIMER(&time_tree[4]);
        
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
#ifdef OPENACC_ENABLED
        #pragma acc exit data delete(sources->x, sources->y, sources->z, sources->q)
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc exit data delete(sources->w)
        }
#endif
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
                total_num_approx += sum_int(let_interaction_list->num_approx, remote_batches->numnodes);
                total_num_direct += sum_int(let_interaction_list->num_direct, remote_batches->numnodes);

                for (int i = 0; i < remote_batches->numnodes; ++i) {
                    for (int j = 0; j < let_interaction_list->num_direct[i]; ++j) {
                        total_num_direct_interact += (long long int) remote_batches->numpar[i]
                            * (long long int) tree->numpar[let_interaction_list->direct_interactions[i][j]];
                    }
                    for (int j = 0; j < let_interaction_list->num_approx[i]; ++j) {
                        total_num_approx_interact += (long long int) remote_batches->numpar[i]
                            * (long long int) run_params->interp_pts_per_cluster;
                    }
                }
            }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            START_TIMER(&time_tree[7]);
#ifdef OPENACC_ENABLED
            #pragma acc enter data copyin(remote_sources->x[0:remote_sources->num], remote_sources->y[0:remote_sources->num], \
                                          remote_sources->z[0:remote_sources->num], remote_sources->q[0:remote_sources->num])
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc enter data create(remote_sources->w[0:remote_sources->num])
            }
#endif
            InteractionCompute_CP(potential, tree, remote_batches, let_interaction_list,
                                  remote_sources, targets, clusters, run_params);
            InteractionLists_Free(&let_interaction_list);
#ifdef OPENACC_ENABLED
            #pragma acc exit data delete(remote_sources->x, remote_sources->y, \
                                         remote_sources->z, remote_sources->q)
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc exit data delete(remote_sources->w)
            }
#endif
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
        // CORRECT AND REORDER POTENTIAL AT TARGETS
        //-------------------------------
        //-------------------------------

        START_TIMER(&time_tree[9]);
        #pragma acc exit data copyout(potential[0:targets->num])
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
#ifdef OPENACC_ENABLED
        #pragma acc exit data delete(targets->x, targets->y, targets->z, \
                                     clusters->x, clusters->y, \
                                     clusters->z, clusters->q)
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc exit data delete(targets->q, clusters->w)
        }
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
#ifdef OPENACC_ENABLED
        #pragma acc enter data copyin(sources->x[0:sources->num], sources->y[0:sources->num], \
                                      sources->z[0:sources->num], sources->q[0:sources->num])
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc enter data copyin(sources->w[0:sources->num])
        }
#endif
        STOP_TIMER(&time_tree[0]);
        
        START_TIMER(&time_tree[1]);
        Batches_Targets_Construct(&batches, targets, run_params);
#ifdef OPENACC_ENABLED
        #pragma acc enter data copyin(targets->x[0:targets->num], targets->y[0:targets->num], \
                                      targets->z[0:targets->num])
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc enter data copyin(targets->q[0:targets->num])
        }
        #pragma acc enter data create(potential[0:targets->num])
#endif
        STOP_TIMER(&time_tree[1]);
        
        START_TIMER(&time_tree[2]);
        Clusters_Sources_Construct(&clusters, sources, tree, run_params);
        STOP_TIMER(&time_tree[2]);
        
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
#ifdef OPENACC_ENABLED
            #pragma acc update self(clusters->x[0:clusters->num], clusters->y[0:clusters->num], \
                                    clusters->z[0:clusters->num], clusters->q[0:clusters->num_charges])
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc update self(clusters->w[0:clusters->num_weights])
            }
#endif
            CommTypesAndTrees_Construct(&comm_types, &let_trees, tree, batches, run_params);
            Particles_Alloc(&let_sources, comm_types->let_sources_length);
            Clusters_Alloc(&let_clusters, comm_types->let_clusters_length, run_params);
                                                    
            CommWindows_Create(&comm_windows, clusters, sources, run_params);
            
            for (int proc_id = 1; proc_id < num_procs; ++proc_id) {
            
                int get_from = (num_procs + rank - proc_id) % num_procs;
                
                CommWindows_Lock(comm_windows, get_from, run_params);
                // This is a non-blocking call!
                CommWindows_GetData(let_clusters, let_sources, comm_types, comm_windows, get_from, run_params);
                CommWindows_Unlock(comm_windows, get_from, run_params);
            }
            
            CommWindows_Free(&comm_windows, run_params);
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
        InteractionCompute_PC(potential, tree, batches, local_interaction_list,
                              sources, targets, clusters, run_params);
        InteractionLists_Free(&local_interaction_list);
#ifdef OPENACC_ENABLED
        #pragma acc exit data delete(sources->x, sources->y, sources->z, sources->q, \
                                     clusters->x, clusters->y, clusters->z, clusters->q)
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc exit data delete(sources->w, clusters->w)
        }
#endif
        STOP_TIMER(&time_tree[5]);

        //~~~~~~~~~~~~~~~~~~~~
        // Remote compute
        //~~~~~~~~~~~~~~~~~~~~
        
        time_tree[6] = 0;
        time_tree[7] = 0;

#ifdef OPENACC_ENABLED
        if (num_procs > 1) {
            START_TIMER(&time1);
            #pragma acc enter data copyin(let_sources->x[0:let_sources->num], let_sources->y[0:let_sources->num], \
                                          let_sources->z[0:let_sources->num], let_sources->q[0:let_sources->num], \
                                          let_clusters->x[0:let_clusters->num], let_clusters->y[0:let_clusters->num], \
                                          let_clusters->z[0:let_clusters->num], let_clusters->q[0:let_clusters->num_charges])
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc enter data copyin(let_sources->w[0:let_sources->num], let_clusters->w[0:let_clusters->num_weights])
            }
            STOP_TIMER(&time1);
            time_tree[6] += time1;
        }
#endif

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

                for (int i = 0; i < batches->numnodes; ++i) {
                    for (int j = 0; j < let_interaction_list->num_direct[i]; ++j) {
                        total_num_direct_interact += (long long int) batches->numpar[i]
                            * (long long int) let_trees[get_from]->numpar[let_interaction_list->direct_interactions[i][j]];
                    }
                    for (int j = 0; j < let_interaction_list->num_approx[i]; ++j) {
                        total_num_approx_interact += (long long int) batches->numpar[i]
                            * (long long int) run_params->interp_pts_per_cluster;
                    }
                }
            }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            START_TIMER(&time1);
            InteractionCompute_PC(potential, let_trees[get_from], batches, let_interaction_list,
                                  let_sources, targets, let_clusters, run_params);
            InteractionLists_Free(&let_interaction_list);
            STOP_TIMER(&time1);
            time_tree[7] += time1;
            
        }

#ifdef OPENACC_ENABLED
        if (num_procs > 1) {
            START_TIMER(&time1);
            #pragma acc exit data delete(let_sources->x, let_sources->y, \
                                         let_sources->z, let_sources->q, \
                                         let_clusters->x, let_clusters->y, \
                                         let_clusters->z, let_clusters->q)
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc exit data delete(let_sources->w, let_clusters->w)
            }
            STOP_TIMER(&time1);
            time_tree[6] += time1;
        }
#endif
            

        //-------------------------------
        //-------------------------------
        // CORRECT AND REORDER
        //-------------------------------
        //-------------------------------
        
        time_tree[8] = 0.0;
        
        START_TIMER(&time_tree[9]);
        #pragma acc exit data copyout(potential[0:targets->num])
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
#ifdef OPENACC_ENABLED
        #pragma acc exit data delete(targets->x, targets->y, targets->z)
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc exit data delete(targets->q)
        }
#endif
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
#ifdef OPENACC_ENABLED
        #pragma acc enter data copyin(sources->x[0:sources->num], sources->y[0:sources->num], \
                                      sources->z[0:sources->num], sources->q[0:sources->num])
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc enter data copyin(sources->w[0:sources->num])
        }
#endif
        STOP_TIMER(&time_tree[0]);

        START_TIMER(&time_tree[1]);
        Tree_Targets_Construct(&target_tree, targets, run_params);
#ifdef OPENACC_ENABLED
        #pragma acc enter data copyin(targets->x[0:targets->num], targets->y[0:targets->num], \
                                      targets->z[0:targets->num])
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc enter data copyin(targets->q[0:targets->num])
        }
        #pragma acc enter data create(potential[0:targets->num])
#endif
        STOP_TIMER(&time_tree[1]);
         
        START_TIMER(&time_tree[2]);
        Clusters_Sources_Construct(&source_clusters, sources, source_tree, run_params);
        Clusters_Targets_Construct(&target_clusters, targets, target_tree, run_params);
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
#ifdef OPENACC_ENABLED
            #pragma acc update self(source_clusters->x[0:source_clusters->num], source_clusters->y[0:source_clusters->num], \
                                    source_clusters->z[0:source_clusters->num], source_clusters->q[0:source_clusters->num_charges])
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc update self(source_clusters->w[0:source_clusters->num_weights])
            }
#endif
            CommTypesAndTrees_Construct(&comm_types, &let_trees,
                                        source_tree, target_tree, run_params);

            Particles_Alloc(&let_sources, comm_types->let_sources_length);
            Clusters_Alloc(&let_clusters, comm_types->let_clusters_length, run_params);
                                                      
            CommWindows_Create(&comm_windows, source_clusters, sources, run_params);
            
            for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

                int get_from = (num_procs + rank - proc_id) % num_procs;
                
                CommWindows_Lock(comm_windows, get_from, run_params);
                //This is a non-blocking call!
                CommWindows_GetData(let_clusters, let_sources, comm_types, comm_windows, get_from, run_params);
                CommWindows_Unlock(comm_windows, get_from, run_params);
            }
            
            CommWindows_Free(&comm_windows, run_params);
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
#ifdef OPENACC_ENABLED
        #pragma acc exit data delete(sources->x, sources->y, sources->z, sources->q, \
                                     source_clusters->x, source_clusters->y, source_clusters->z, source_clusters->q)
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc exit data delete(sources->w, source_clusters->w)
        }
#endif
        STOP_TIMER(&time_tree[5]);
        
        //~~~~~~~~~~~~~~~~~~~~
        // Remote compute
        //~~~~~~~~~~~~~~~~~~~~
        
        time_tree[6] = 0;
        time_tree[7] = 0;

#ifdef OPENACC_ENABLED
        if (num_procs > 1) {
            START_TIMER(&time1);
            #pragma acc enter data copyin(let_sources->x[0:let_sources->num], let_sources->y[0:let_sources->num], \
                                          let_sources->z[0:let_sources->num], let_sources->q[0:let_sources->num], \
                                          let_clusters->x[0:let_clusters->num], let_clusters->y[0:let_clusters->num], \
                                          let_clusters->z[0:let_clusters->num], let_clusters->q[0:let_clusters->num_charges])
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc enter data copyin(let_sources->w[0:let_sources->num], let_clusters->w[0:let_clusters->num_weights])
            }
            STOP_TIMER(&time1);
            time_tree[6] += time1;
        }
#endif
            
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
                
                total_num_source_approx += sum_int(let_interaction_list->num_cc_source_approx,
                                                   target_tree->numnodes);
                total_num_target_approx += sum_int(let_interaction_list->num_cc_target_approx,
                                                   target_tree->numnodes);


                total_num_approx_interact += (long long int) sum_int(let_interaction_list->num_approx, target_tree->numnodes)
                                           * (long long int) run_params->interp_pts_per_cluster
                                           * (long long int) run_params->interp_pts_per_cluster;

                for (int i = 0; i < target_tree->numnodes; ++i) {
                    for (int j = 0; j < let_interaction_list->num_direct[i]; ++j) {
                        total_num_direct_interact += (long long int) target_tree->numpar[i]
                            * (long long int) let_trees[get_from]->numpar[let_interaction_list->direct_interactions[i][j]];
                    }
                    for (int j = 0; j < let_interaction_list->num_cc_source_approx[i]; ++j) {
                        total_num_source_approx_interact += (long long int) target_tree->numpar[i]
                            * (long long int) run_params->interp_pts_per_cluster;
                    }
                    for (int j = 0; j < let_interaction_list->num_cc_target_approx[i]; ++j) {
                        total_num_target_approx_interact += (long long int) run_params->interp_pts_per_cluster
                            * (long long int) let_trees[get_from]->numpar[let_interaction_list->cc_target_approx_interactions[i][j]];
                    }
                }
            }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~


            START_TIMER(&time1);
            InteractionCompute_CC(potential, let_trees[get_from], target_tree, let_interaction_list,
                                  let_sources, targets, let_clusters, target_clusters, run_params);
            InteractionLists_Free(&let_interaction_list);
            STOP_TIMER(&time1);
            time_tree[7] += time1;
        }

#ifdef OPENACC_ENABLED
        if (num_procs > 1) {
            START_TIMER(&time1);
            #pragma acc exit data delete(let_sources->x, let_sources->y, \
                                         let_sources->z, let_sources->q, \
                                         let_clusters->x, let_clusters->y, \
                                         let_clusters->z, let_clusters->q)
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc exit data delete(let_sources->w, let_clusters->w)
            }
            STOP_TIMER(&time1);
            time_tree[6] += time1;
        }
#endif
            
            
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
        #pragma acc exit data copyout(potential[0:targets->num])
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
#ifdef OPENACC_ENABLED
        #pragma acc exit data delete(targets->x, targets->y, targets->z, \
                                     target_clusters->x, target_clusters->y, \
                                     target_clusters->z, target_clusters->q)
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc exit data delete(targets->q, target_clusters->w)
        }
#endif
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


    if (sources_w_dummy) free_vector(sources->w);
    if (targets_q_dummy) free_vector(targets->q);


//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    if (run_params->verbosity > 0) {
       
        int global_num_inter,  max_num_inter,  min_num_inter;
        int global_num_direct, max_num_direct, min_num_direct;
        int global_num_approx, max_num_approx, min_num_approx;
        
        int global_num_source_approx, max_num_source_approx, min_num_source_approx;
        int global_num_target_approx, max_num_target_approx, min_num_target_approx;

        total_num_inter = total_num_direct + total_num_approx
                        + total_num_source_approx + total_num_target_approx;
                        
        MPI_Reduce(&total_num_inter,   &global_num_inter, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_inter,      &max_num_inter, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_inter,      &min_num_inter, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&total_num_direct, &global_num_direct, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_direct,    &max_num_direct, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_direct,    &min_num_direct, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&total_num_approx, &global_num_approx, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_approx,    &max_num_approx, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_approx,    &min_num_approx, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        // These types of interactions only occur for CC
        if (run_params->compute_type == CLUSTER_CLUSTER) {
            MPI_Reduce(&total_num_source_approx, &global_num_source_approx, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&total_num_source_approx,    &max_num_source_approx, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&total_num_source_approx,    &min_num_source_approx, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
            MPI_Reduce(&total_num_target_approx, &global_num_target_approx, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&total_num_target_approx,    &max_num_target_approx, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&total_num_target_approx,    &min_num_target_approx, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        }
        
        if (rank == 0) {
            printf("[BaryTree]\n");
            printf("[BaryTree] Interaction information: \n");
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
            
            // These types of interactions only occur for CC
            if (run_params->compute_type == CLUSTER_CLUSTER) {
                printf("[BaryTree] Cumulative source approx inter across all ranks: %d\n", global_num_source_approx);
                printf("[BaryTree]    Maximum source approx inter across all ranks: %d\n", max_num_source_approx);
                printf("[BaryTree]    Minimum source approx inter across all ranks: %d\n", min_num_source_approx);
                printf("[BaryTree]                                           Ratio: %f\n",
                       (double)max_num_source_approx / (double)min_num_source_approx);
                printf("[BaryTree]\n");
                printf("[BaryTree] Cumulative target approx inter across all ranks: %d\n", global_num_target_approx);
                printf("[BaryTree]    Maximum target approx inter across all ranks: %d\n", max_num_target_approx);
                printf("[BaryTree]    Minimum target approx inter across all ranks: %d\n", min_num_target_approx);
                printf("[BaryTree]                                           Ratio: %f\n",
                       (double)max_num_target_approx / (double)min_num_target_approx);
                printf("[BaryTree]\n");
            }
        }


        /* For the pointwise interactions */

        long long int global_num_interact,  max_num_interact,  min_num_interact;
        long long int global_num_direct_interact, max_num_direct_interact, min_num_direct_interact;
        long long int global_num_approx_interact, max_num_approx_interact, min_num_approx_interact;

        long long int global_num_source_approx_interact, max_num_source_approx_interact, min_num_source_approx_interact;
        long long int global_num_target_approx_interact, max_num_target_approx_interact, min_num_target_approx_interact;

        total_num_interact = total_num_direct_interact + total_num_approx_interact
                           + total_num_source_approx_interact + total_num_target_approx_interact;
                        
        MPI_Reduce(&total_num_interact,   &global_num_interact, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_interact,      &max_num_interact, 1, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_interact,      &min_num_interact, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&total_num_direct_interact, &global_num_direct_interact, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_direct_interact,    &max_num_direct_interact, 1, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_direct_interact,    &min_num_direct_interact, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&total_num_approx_interact, &global_num_approx_interact, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_approx_interact,    &max_num_approx_interact, 1, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&total_num_approx_interact,    &min_num_approx_interact, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        // These types of interactions only occur for CC
        if (run_params->compute_type == CLUSTER_CLUSTER) {
            MPI_Reduce(&total_num_source_approx_interact, &global_num_source_approx_interact, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&total_num_source_approx_interact,    &max_num_source_approx_interact, 1, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&total_num_source_approx_interact,    &min_num_source_approx_interact, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
            MPI_Reduce(&total_num_target_approx_interact, &global_num_target_approx_interact, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            MPI_Reduce(&total_num_target_approx_interact,    &max_num_target_approx_interact, 1, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
            MPI_Reduce(&total_num_target_approx_interact,    &min_num_target_approx_interact, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        }

        if (rank == 0) {
            printf("[BaryTree]\n");
            printf("[BaryTree]               Cumulative pointwise interactions across all ranks: %lld\n", global_num_interact);
            printf("[BaryTree]                  Maximum pointwise interactions across all ranks: %lld\n", max_num_interact);
            printf("[BaryTree]                  Minimum pointwise interactions across all ranks: %lld\n", min_num_interact);
            printf("[BaryTree]\n"); 

            printf("[BaryTree]        Cumulative direct pointwise interactions across all ranks: %lld\n", global_num_direct_interact);
            printf("[BaryTree]           Maximum direct pointwise interactions across all ranks: %lld\n", max_num_direct_interact);
            printf("[BaryTree]           Minimum direct pointwise interactions across all ranks: %lld\n", min_num_direct_interact);
            printf("[BaryTree]\n");

            printf("[BaryTree]        Cumulative approx pointwise interactions across all ranks: %lld\n", global_num_approx_interact);
            printf("[BaryTree]           Maximum approx pointwise interactions across all ranks: %lld\n", max_num_approx_interact);
            printf("[BaryTree]           Minimum approx pointwise interactions across all ranks: %lld\n", min_num_approx_interact);
            printf("[BaryTree]\n"); 

            // These types of interactions only occur for CC
            if (run_params->compute_type == CLUSTER_CLUSTER) {
                printf("[BaryTree] Cumulative source approx pointwise interactions across all ranks: %lld\n", global_num_source_approx_interact);
                printf("[BaryTree]    Maximum source approx pointwise interactions across all ranks: %lld\n", max_num_source_approx_interact);
                printf("[BaryTree]    Minimum source approx pointwise interactions across all ranks: %lld\n", min_num_source_approx_interact);
                printf("[BaryTree]\n");

                printf("[BaryTree] Cumulative target approx pointwise interactions across all ranks: %lld\n", global_num_target_approx_interact);
                printf("[BaryTree]    Maximum source approx pointwise interactions across all ranks: %lld\n", max_num_target_approx_interact);
                printf("[BaryTree]    Minimum source approx pointwise interactions across all ranks: %lld\n", min_num_target_approx_interact);
                printf("[BaryTree]\n");
            }
            
            printf("[BaryTree] BaryTree has finished.\n");
            printf("[BaryTree]\n");
        }
    }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    
    return;
}
