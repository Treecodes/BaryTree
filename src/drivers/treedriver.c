#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <limits.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../utilities/timers.h"
#include "../utilities/enums.h"
#include "../utilities/advanced_timings.h"

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
    Particles_Validate(sources, targets, run_params);
    Particles_ConstructOrder(sources);
    Particles_ConstructOrder(targets);

    
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
    if (run_params->verbosity > 0 && rank == 0) {
        printf("[BaryTree]\n");
        printf("[BaryTree] BaryTree has started.\n");
        RunParams_Print(run_params);
    }
    if (run_params->verbosity > 1 ) {


        int M_min_g, M_max_g, M_avg_g, N_min_g, N_max_g, N_avg_g;

        int M_max = targets->num;
        int M_min = targets->num;
        int M_avg = targets->num;

        int N_max = sources->num;
        int N_min = sources->num;
        int N_avg = sources->num;

        MPI_Reduce(&M_max, &M_max_g, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&M_min, &M_min_g, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&M_avg, &M_avg_g, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        M_avg_g = M_avg_g / num_procs;

        MPI_Reduce(&N_max, &N_max_g, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&N_min, &N_min_g, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        MPI_Reduce(&N_avg, &N_avg_g, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        N_avg_g = N_avg_g / num_procs;

        if ( rank == 0) {
            printf("[BaryTree] %d target and %d source particles distributed over %d MPI ranks.\n", M_avg_g*num_procs, N_avg_g*num_procs, num_procs);
            printf("[BaryTree] min, max, avg targets per rank: (%d,%d,%d).\n", M_min_g,M_max_g,M_avg_g);
            printf("[BaryTree] min, max, avg sources per rank: (%d,%d,%d).\n", N_min_g,N_max_g,N_avg_g);
        }
    }
    
    double time1;
    long long num_pp = 0;
    long long num_cc = 0;
    long long num_pc = 0;
    long long num_cp = 0;

    long long num_pp_ptwise = 0;
    long long num_cc_ptwise = 0;
    long long num_pc_ptwise = 0;
    long long num_cp_ptwise = 0;
    
    long long num_cc_replaced = 0;
    long long num_pc_replaced = 0;
    long long num_cp_replaced = 0;

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
        #pragma acc kernels present(potential)
        {
            for (int i=0;i<targets->num;i++){
                potential[i]=0.0;
            }
        }
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
        if (run_params->verbosity > 2) {
            num_cp += sum_int(local_interaction_list->num_cp, batches->numnodes);
            num_pp += sum_int(local_interaction_list->num_pp, batches->numnodes);

            for (int i = 0; i < batches->numnodes; ++i) {
                for (int j = 0; j < local_interaction_list->num_pp[i]; ++j) {
                    num_pp_ptwise += (long long) batches->numpar[i]
                        * (long long) tree->numpar[local_interaction_list->pp_interactions[i][j]];
                }
                for (int j = 0; j < local_interaction_list->num_cp[i]; ++j) {
                    num_cp_ptwise += (long long) batches->numpar[i]
                        * (long long) run_params->interp_pts_per_cluster;

                    num_cp_replaced += (long long) batches->numpar[i]
                        * (long long) tree->numpar[local_interaction_list->cp_interactions[i][j]];
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
            if (run_params->verbosity > 2) {
                num_cp += sum_int(let_interaction_list->num_cp, remote_batches->numnodes);
                num_pp += sum_int(let_interaction_list->num_pp, remote_batches->numnodes);

                for (int i = 0; i < remote_batches->numnodes; ++i) {
                    for (int j = 0; j < let_interaction_list->num_pp[i]; ++j) {
                        num_pp_ptwise += (long long) remote_batches->numpar[i]
                            * (long long) tree->numpar[let_interaction_list->pp_interactions[i][j]];
                    }
                    for (int j = 0; j < let_interaction_list->num_cp[i]; ++j) {
                        num_cp_ptwise += (long long) remote_batches->numpar[i]
                            * (long long) run_params->interp_pts_per_cluster;
                            
                        num_cp_replaced += (long long) remote_batches->numpar[i]
                            * (long long) tree->numpar[let_interaction_list->cp_interactions[i][j]];
                    }
                }
            }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

            START_TIMER(&time_tree[7]);
#ifdef OPENACC_ENABLED
            #pragma acc enter data copyin(remote_sources->x[0:remote_sources->num], remote_sources->y[0:remote_sources->num], \
                                          remote_sources->z[0:remote_sources->num], remote_sources->q[0:remote_sources->num])
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc enter data copyin(remote_sources->w[0:remote_sources->num])


            }
#endif
            InteractionCompute_CP(potential, tree, remote_batches, let_interaction_list,
                                  remote_sources, targets, clusters, run_params);
            InteractionLists_Free(&let_interaction_list);
#ifdef OPENACC_ENABLED
            #pragma acc exit data delete(remote_sources->x[0:remote_sources->num], remote_sources->y[0:remote_sources->num], \
                                         remote_sources->z[0:remote_sources->num], remote_sources->q[0:remote_sources->num])
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc exit data delete(remote_sources->w[0:remote_sources->num])
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
        #pragma acc exit data delete(targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num], \
                                     clusters->x[0:clusters->num], clusters->y[0:clusters->num], \
                                     clusters->z[0:clusters->num], clusters->q[0:clusters->num_charges])
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc exit data delete(targets->q[0:targets->num], clusters->w[0:clusters->num_weights])
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

        #pragma acc kernels present(potential)
        {
            for (int i=0;i<targets->num;i++){
                potential[i]=0.0;
            }
        }
#endif
        STOP_TIMER(&time_tree[1]);
        
        START_TIMER(&time_tree[2]);
        Clusters_Sources_Construct(&clusters, sources, tree, run_params);
        STOP_TIMER(&time_tree[2]);
        
        //~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
        if (run_params->verbosity > 1) {
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
        if (run_params->verbosity > 2) {
            num_pc += sum_int(local_interaction_list->num_pc, batches->numnodes);
            num_pp += sum_int(local_interaction_list->num_pp, batches->numnodes);

            for (int i = 0; i < batches->numnodes; ++i) {
                for (int j = 0; j < local_interaction_list->num_pp[i]; ++j) {
                    num_pp_ptwise += (long long) batches->numpar[i]
                        * (long long) tree->numpar[local_interaction_list->pp_interactions[i][j]];
                }
                
                for (int j = 0; j < local_interaction_list->num_pc[i]; ++j) {
                    num_pc_ptwise += (long long) batches->numpar[i]
                        * (long long) run_params->interp_pts_per_cluster;
                        
                    num_pc_replaced += (long long) batches->numpar[i]
                        * (long long) tree->numpar[local_interaction_list->pc_interactions[i][j]];
                }
            }
        }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
  
        START_TIMER(&time_tree[5]);
        InteractionCompute_PC(potential, tree, batches, local_interaction_list,
                              sources, targets, clusters, run_params);
        InteractionLists_Free(&local_interaction_list);
#ifdef OPENACC_ENABLED
        #pragma acc exit data delete(sources->x[0:sources->num], sources->y[0:sources->num], sources->z[0:sources->num], sources->q[0:sources->num], \
                                     clusters->x[0:clusters->num], clusters->y[0:clusters->num], clusters->z[0:clusters->num], \
                                     clusters->q[0:clusters->num_charges])
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc exit data delete(sources->w[0:sources->num], clusters->w[0:clusters->num_weights])
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
            if (run_params->verbosity > 2) {
                num_pc += sum_int(let_interaction_list->num_pc, batches->numnodes);
                num_pp += sum_int(let_interaction_list->num_pp, batches->numnodes);

                for (int i = 0; i < batches->numnodes; ++i) {
                    for (int j = 0; j < let_interaction_list->num_pp[i]; ++j) {
                        num_pp_ptwise += (long long) batches->numpar[i]
                            * (long long) let_trees[get_from]->numpar[let_interaction_list->pp_interactions[i][j]];
                    }
                    
                    for (int j = 0; j < let_interaction_list->num_pc[i]; ++j) {
                        num_pc_ptwise += (long long) batches->numpar[i]
                            * (long long) run_params->interp_pts_per_cluster;
                            
                        num_pc_replaced += (long long) batches->numpar[i]
                            * (long long) let_trees[get_from]->numpar[let_interaction_list->pc_interactions[i][j]];
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
            #pragma acc exit data delete(let_sources->x[0:let_sources->num], let_sources->y[0:let_sources->num], \
                                         let_sources->z[0:let_sources->num], let_sources->q[0:let_sources->num], \
                                         let_clusters->x[0:let_clusters->num], let_clusters->y[0:let_clusters->num], \
                                         let_clusters->z[0:let_clusters->num], let_clusters->q[0:let_clusters->num_charges])
            if (run_params->singularity == SUBTRACTION) {
                #pragma acc exit data delete(let_sources->w[0:let_sources->num], let_clusters->w[0:let_clusters->num_weights])
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
        #pragma acc exit data delete(targets->x[0:targets->num], targets->y[0:targets->num], targets->z[0:targets->num])
        if (run_params->singularity == SUBTRACTION) {
            #pragma acc exit data delete(targets->q[0:targets->num])
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
        #pragma acc kernels present(potential)
        {
            for (int i=0;i<targets->num;i++){
                potential[i]=0.0;
            }
        }
#endif
        STOP_TIMER(&time_tree[1]);
         
        START_TIMER(&time_tree[2]);
        Clusters_Sources_Construct(&source_clusters, sources, source_tree, run_params);
        Clusters_Targets_Construct(&target_clusters, targets, target_tree, run_params);
        STOP_TIMER(&time_tree[2]);
        
        //~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
                if (run_params->verbosity > 1) {
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
        if (run_params->verbosity > 2) {
            num_cc += sum_int(local_interaction_list->num_cc, target_tree->numnodes);
            num_pp += sum_int(local_interaction_list->num_pp, target_tree->numnodes);
            num_pc += sum_int(local_interaction_list->num_pc, target_tree->numnodes);
            num_cp += sum_int(local_interaction_list->num_cp, target_tree->numnodes);

            num_cc_ptwise += (long long) sum_int(local_interaction_list->num_cc, target_tree->numnodes)
                                   * (long long) run_params->interp_pts_per_cluster
                                   * (long long) run_params->interp_pts_per_cluster;
  
            for (int i = 0; i < target_tree->numnodes; ++i) {
                for (int j = 0; j < local_interaction_list->num_pp[i]; ++j) {
                    num_pp_ptwise += (long long) target_tree->numpar[i]
                        * (long long) source_tree->numpar[local_interaction_list->pp_interactions[i][j]];
                }
                
                for (int j = 0; j < local_interaction_list->num_cc[i]; ++j) {
                    num_cc_replaced += (long long) target_tree->numpar[i]
                        * (long long) source_tree->numpar[local_interaction_list->cc_interactions[i][j]];
                }
                
                for (int j = 0; j < local_interaction_list->num_pc[i]; ++j) {
                    num_pc_ptwise += (long long) target_tree->numpar[i]
                        * (long long) run_params->interp_pts_per_cluster;
                        
                    num_pc_replaced += (long long) target_tree->numpar[i]
                        * (long long) source_tree->numpar[local_interaction_list->pc_interactions[i][j]];
                }
                
                for (int j = 0; j < local_interaction_list->num_cp[i]; ++j) {
                    num_cp_ptwise += (long long) run_params->interp_pts_per_cluster
                        * (long long) source_tree->numpar[local_interaction_list->cp_interactions[i][j]];
                        
                    num_cp_replaced += (long long) target_tree->numpar[i]
                        * (long long) source_tree->numpar[local_interaction_list->cp_interactions[i][j]];
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
            if (run_params->verbosity > 2) {
                num_pp += sum_int(let_interaction_list->num_pp, target_tree->numnodes);
                num_cc += sum_int(let_interaction_list->num_cc, target_tree->numnodes);
                num_pc += sum_int(let_interaction_list->num_pc, target_tree->numnodes);
                num_cp += sum_int(let_interaction_list->num_cp, target_tree->numnodes);

                num_cc_ptwise += (long long) sum_int(let_interaction_list->num_cc, target_tree->numnodes)
                                       * (long long) run_params->interp_pts_per_cluster
                                       * (long long) run_params->interp_pts_per_cluster;

                for (int i = 0; i < target_tree->numnodes; ++i) {
                    for (int j = 0; j < let_interaction_list->num_pp[i]; ++j) {
                        num_pp_ptwise += (long long) target_tree->numpar[i]
                            * (long long) let_trees[get_from]->numpar[let_interaction_list->pp_interactions[i][j]];
                    }
                    
                    for (int j = 0; j < let_interaction_list->num_cc[i]; ++j) {
                        num_cc_replaced += (long long) target_tree->numpar[i]
                            * (long long) let_trees[get_from]->numpar[let_interaction_list->cc_interactions[i][j]];
                    }
                    
                    for (int j = 0; j < let_interaction_list->num_pc[i]; ++j) {
                        num_pc_ptwise += (long long) target_tree->numpar[i]
                            * (long long) run_params->interp_pts_per_cluster;
                            
                        num_pc_replaced += (long long) target_tree->numpar[i]
                            * (long long) let_trees[get_from]->numpar[let_interaction_list->pc_interactions[i][j]];
                    }
                    
                    for (int j = 0; j < let_interaction_list->num_cp[i]; ++j) {
                        num_cp_ptwise += (long long) run_params->interp_pts_per_cluster
                            * (long long) let_trees[get_from]->numpar[let_interaction_list->cp_interactions[i][j]];
                            
                        num_cp_replaced += (long long) target_tree->numpar[i]
                            * (long long) let_trees[get_from]->numpar[let_interaction_list->cp_interactions[i][j]];
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


    if (sources->sources_w_dummy) free_vector(sources->w);
    if (targets->targets_q_dummy) free_vector(targets->q);


//~ ~ ~ D I A G N O S T I C S ~ ~ ~ S T A R T ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    double total_time[1], total_time_glob[1];
    if (run_params->verbosity > 0) {

        /* Total treedriver time */
        total_time[0] = time_tree[0] + time_tree[1] + time_tree[2] + time_tree[3] + time_tree[4] +
                            time_tree[5] + time_tree[6] + time_tree[7] + time_tree[8] + time_tree[9] +
                            time_tree[10];

        if (rank==0) {
            printf("[BaryTree] Total BaryTree time: %1.3f seconds.\n", total_time[0]);
        }
    }

    if (run_params->verbosity > 2) {
       
        int global_num_all = 0,  max_num_all = 0,  min_num_all = 0;
        int global_num_pp = 0, max_num_pp = 0, min_num_pp = 0;
        int global_num_cc = 0, max_num_cc = 0, min_num_cc = 0;
        int global_num_pc = 0, max_num_pc = 0, min_num_pc = 0;
        int global_num_cp = 0, max_num_cp = 0, min_num_cp = 0;

        long long num_all = num_pp + num_cc + num_pc + num_cp;
                        
        MPI_Reduce(&num_all, &global_num_all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_all,    &max_num_all, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_all,    &min_num_all, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&num_pp,   &global_num_pp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_pp,      &max_num_pp, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_pp,      &min_num_pp, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&num_cc,   &global_num_cc, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_cc,      &max_num_cc, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_cc,      &min_num_cc, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&num_pc,   &global_num_pc, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_pc,      &max_num_pc, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_pc,      &min_num_pc, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    
        MPI_Reduce(&num_cp,   &global_num_cp, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_cp,      &max_num_cp, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_cp,      &min_num_cp, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        if (rank == 0) {
            printf("[BaryTree]\n");
            printf("[BaryTree] Clusterwise interaction information:\n");
            printf("[BaryTree]\n");
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] |       |             Cumulative       Minimum per rank       Maximum per rank |\n");
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] | Total | %22d %22d %22d |\n", global_num_all, min_num_all, max_num_all);
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] | PP    | %22d %22d %22d |\n", global_num_pp, min_num_pp, max_num_pp);
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] | PC    | %22d %22d %22d |\n", global_num_pc, min_num_pc, max_num_pc);
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] | CP    | %22d %22d %22d |\n", global_num_cp, min_num_cp, max_num_cp);
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] | CC    | %22d %22d %22d |\n", global_num_cc, min_num_cc, max_num_cc);
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree]\n");
        }

        /* For the pointwise interactions */

        long long global_num_all_ptwise = 0,  max_num_all_ptwise = 0,  min_num_all_ptwise = 0;
        long long global_num_pp_ptwise = 0, max_num_pp_ptwise = 0, min_num_pp_ptwise = 0;
        long long global_num_cc_ptwise = 0, max_num_cc_ptwise = 0, min_num_cc_ptwise = 0;
        long long global_num_pc_ptwise = 0, max_num_pc_ptwise = 0, min_num_pc_ptwise = 0;
        long long global_num_cp_ptwise = 0, max_num_cp_ptwise = 0, min_num_cp_ptwise = 0;

        long long num_all_ptwise = num_pp_ptwise + num_cc_ptwise + num_pc_ptwise + num_cp_ptwise;
                        
        MPI_Reduce(&num_all_ptwise,   &global_num_all_ptwise, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_all_ptwise,      &max_num_all_ptwise, 1, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_all_ptwise,      &min_num_all_ptwise, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&num_pp_ptwise, &global_num_pp_ptwise, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_pp_ptwise,    &max_num_pp_ptwise, 1, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_pp_ptwise,    &min_num_pp_ptwise, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&num_cc_ptwise, &global_num_cc_ptwise, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_cc_ptwise,    &max_num_cc_ptwise, 1, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_cc_ptwise,    &min_num_cc_ptwise, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        
        MPI_Reduce(&num_pc_ptwise, &global_num_pc_ptwise, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_pc_ptwise,    &max_num_pc_ptwise, 1, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_pc_ptwise,    &min_num_pc_ptwise, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    
        MPI_Reduce(&num_cp_ptwise, &global_num_cp_ptwise, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_cp_ptwise,    &max_num_cp_ptwise, 1, MPI_LONG_LONG_INT, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_cp_ptwise,    &min_num_cp_ptwise, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);

        if (rank == 0) {
            printf("[BaryTree]\n");
            printf("[BaryTree] Pointwise interaction information:\n");
            printf("[BaryTree]\n");
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] |       |             Cumulative       Minimum per rank       Maximum per rank |\n");
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] | Total | %22lld %22lld %22lld |\n",
                   global_num_all_ptwise, min_num_all_ptwise, max_num_all_ptwise);
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] | PP    | %22lld %22lld %22lld |\n",
                   global_num_pp_ptwise, min_num_pp_ptwise, max_num_pp_ptwise);
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] | PC    | %22lld %22lld %22lld |\n",
                   global_num_pc_ptwise, min_num_pc_ptwise, max_num_pc_ptwise);
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] | CP    | %22lld %22lld %22lld |\n",
                   global_num_cp_ptwise, min_num_cp_ptwise, max_num_cp_ptwise);
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree] | CC    | %22lld %22lld %22lld |\n",
                   global_num_cc_ptwise, min_num_cc_ptwise, max_num_cc_ptwise);
            printf("[BaryTree] |------------------------------------------------------------------------------|\n");
            printf("[BaryTree]\n");
        }
        
        long long global_num_cc_replaced = 0;
        long long global_num_pc_replaced = 0;
        long long global_num_cp_replaced = 0;
        
        MPI_Reduce(&num_cc_replaced, &global_num_cc_replaced, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_pc_replaced, &global_num_pc_replaced, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&num_cp_replaced, &global_num_cp_replaced, 1, MPI_LONG_LONG_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
        long long global_num_pure_pp = global_num_pp_ptwise   + global_num_cc_replaced
                                     + global_num_pc_replaced + global_num_cp_replaced;
                           
        if (rank == 0) {
            printf("[BaryTree]\n");
            printf("[BaryTree] Pointwise interactions accounted for:\n");
            printf("[BaryTree]\n");
            printf("[BaryTree] |-------------------------------------------------------|\n");
            printf("[BaryTree] |       |                 Number             Percentage |\n");
            printf("[BaryTree] |-------------------------------------------------------|\n");
            printf("[BaryTree] | Total | %22lld            100.000000%% |\n",
                   global_num_pure_pp);
            printf("[BaryTree] |-------------------------------------------------------|\n");
            printf("[BaryTree] | PP    | %22lld %21.6f%% |\n",
                   global_num_pp_ptwise,   100. * (double)global_num_pp_ptwise / (double)global_num_pure_pp);
            printf("[BaryTree] |-------------------------------------------------------|\n");
            printf("[BaryTree] | PC    | %22lld %21.6f%% |\n",
                   global_num_pc_replaced, 100. * (double)global_num_pc_replaced / (double)global_num_pure_pp);
            printf("[BaryTree] |-------------------------------------------------------|\n");
            printf("[BaryTree] | CP    | %22lld %21.6f%% |\n",
                   global_num_cp_replaced, 100. * (double)global_num_cp_replaced / (double)global_num_pure_pp);
            printf("[BaryTree] |-------------------------------------------------------|\n");
            printf("[BaryTree] | CC    | %22lld %21.6f%% |\n",
                   global_num_cc_replaced, 100. * (double)global_num_cc_replaced / (double)global_num_pure_pp);
            printf("[BaryTree] |-------------------------------------------------------|\n");
            printf("[BaryTree]\n");
            
        }

        /* variables for date-time calculation */
        double time_tree_glob[3][13];

        Timing_Calculate(time_tree_glob, time_tree, total_time_glob, total_time);
        Timing_Print(time_tree_glob, total_time_glob, run_params);
    }
//~ ~ ~ D I A G N O S T I C S ~ ~ ~ E N D ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

    if (rank == 0) {
            printf("[BaryTree] BaryTree has finished.\n");
            printf("[BaryTree]\n");
    }

    return;
}
