#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <limits.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../utilities/enums.h"

#include "../globvars.h"

#include "../tree/struct_nodes.h"
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

#include "../run_params/struct_run_params.h"
#include "../run_params/run_params.h"

#include "../interaction_lists/interaction_lists.h"
#include "../interaction_compute/interaction_compute.h"

#include "treedriver.h"


void treedriver(struct particles *sources, struct particles *targets, struct RunParams *run_params,
                double *potential_array, double *time_tree)
{
    int rank, num_procs, ierr;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    RunParams_Validate(run_params);
    if (run_params->verbosity > 0) printf("Set rank %i and num_procs %i.\n", rank, num_procs);

    double time1;
    

    int totalNumberDirect = 0;
    int totalNumberApprox = 0;
    int totalNumberInteractions = 0;
    int cumulativeNumberInteractions = 0;
    int maxNumberInteractions = 0;
    int minNumberInteractions = 0;



    //--------------------------------------
    //--------------------------------------
    //CLUSTER PARTICLE
    //--------------------------------------
    //--------------------------------------
    


    /* call setup to allocate arrays for Taylor expansions and setup global vars */
    if (run_params->compute_type == CLUSTER_PARTICLE) {
    
        struct tnode *troot = NULL;
        struct tnode_array *tree_array = NULL;
        double xyzminmax[6];
        int numnodes = 0, numleaves = 0;

        struct tnode_array *batches = NULL;
        double batch_lim[6];

        struct clusters *clusters = NULL;
        

        time1 = MPI_Wtime();

        Tree_Setup(sources, targets, run_params->interp_order, xyzminmax);
        Tree_CP_Create(&troot, targets, 1, targets->num,
                       run_params->max_per_target_leaf, xyzminmax, 0, &numnodes, &numleaves);
        Tree_SetIndex(troot, 0);
        Tree_AllocArray(&tree_array, numnodes);
        Tree_CreateArray(troot, tree_array);


        time_tree[0] = MPI_Wtime() - time1; //time_maketreearray
        


        time1 = MPI_Wtime();

        Batches_Alloc(&batches, batch_lim, sources, run_params->max_per_source_leaf);
        Batches_CreateSourceBatches(batches, sources, 1, sources->num, run_params->max_per_source_leaf, batch_lim);


        time_tree[1] = MPI_Wtime() - time1; //time_createbatch
        

        time1 = MPI_Wtime();

        Clusters_CP_Setup(&clusters, run_params->interp_order, tree_array,
                          run_params->approximation, run_params->singularity);


        time_tree[2] = MPI_Wtime() - time1; //time_fillclusters
        
        
        //-------------------
        //BEGIN COMPUTE PHASE
        //-------------------
        
        
        MPI_Barrier(MPI_COMM_WORLD);

        time1 = MPI_Wtime();

        int num_batches_on_proc[num_procs];
        MPI_Allgather(&(batches->numnodes), 1, MPI_INT, num_batches_on_proc, 1, MPI_INT, MPI_COMM_WORLD);
        
        int num_sources_on_proc[num_procs];
        MPI_Allgather(&(sources->num), 1, MPI_INT, num_sources_on_proc, 1, MPI_INT, MPI_COMM_WORLD);
        
        int num_remote_batches = sum_int(num_batches_on_proc, num_procs) - num_batches_on_proc[rank];
        int num_remote_sources = sum_int(num_sources_on_proc, num_procs) - num_sources_on_proc[rank];
        
        int previous_remote_sources_length[num_procs];

        struct tnode_array *remote_batches = NULL;
        struct particles *remote_sources = NULL;
        
        Batches_AllocArray(&remote_batches, num_remote_batches);


        MPI_Win win_x_mid, win_y_mid, win_z_mid, win_radius, win_numpar, win_ibeg, win_iend;
        MPI_Win win_sources_x, win_sources_y, win_sources_z, win_sources_q, win_sources_w;
        
        MPI_Win_create(batches->x_mid,  batches->numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_x_mid);
        MPI_Win_create(batches->y_mid,  batches->numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_y_mid);
        MPI_Win_create(batches->z_mid,  batches->numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_z_mid);
        MPI_Win_create(batches->radius, batches->numnodes*sizeof(double), sizeof(double),  MPI_INFO_NULL, MPI_COMM_WORLD, &win_radius);
        MPI_Win_create(batches->numpar, batches->numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_numpar);
        MPI_Win_create(batches->ibeg,   batches->numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_ibeg);
        MPI_Win_create(batches->iend,   batches->numnodes*sizeof(int),    sizeof(int),     MPI_INFO_NULL, MPI_COMM_WORLD, &win_iend);


        MPI_Win_create(sources->x, sources->num*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_x);
        MPI_Win_create(sources->y, sources->num*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_y);
        MPI_Win_create(sources->z, sources->num*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_z);
        MPI_Win_create(sources->q, sources->num*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_q);
        MPI_Win_create(sources->w, sources->num*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_w);

        // Perform MPI round robin, filling LET with remote data
        
        int previous_remote_batches_counter = 0;
        int previous_remote_sources_counter = 0;

        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

            int get_from = (num_procs+rank-proc_id) % num_procs;

            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_x_mid);
            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_y_mid);
            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_z_mid);
            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_radius);
            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_numpar);
            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_ibeg);
            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_iend);

            MPI_Get(&(remote_batches->x_mid[previous_remote_batches_counter]),  num_batches_on_proc[get_from], MPI_DOUBLE,
                    get_from, 0, num_batches_on_proc[get_from], MPI_DOUBLE, win_x_mid);
                    
            MPI_Get(&(remote_batches->y_mid[previous_remote_batches_counter]),  num_batches_on_proc[get_from], MPI_DOUBLE,
                    get_from, 0, num_batches_on_proc[get_from], MPI_DOUBLE, win_y_mid);
                    
            MPI_Get(&(remote_batches->z_mid[previous_remote_batches_counter]),  num_batches_on_proc[get_from], MPI_DOUBLE,
                    get_from, 0, num_batches_on_proc[get_from], MPI_DOUBLE, win_z_mid);
                    
            MPI_Get(&(remote_batches->radius[previous_remote_batches_counter]), num_batches_on_proc[get_from], MPI_DOUBLE,
                    get_from, 0, num_batches_on_proc[get_from], MPI_DOUBLE, win_radius);
                    
            MPI_Get(&(remote_batches->numpar[previous_remote_batches_counter]), num_batches_on_proc[get_from], MPI_INT,
                    get_from, 0, num_batches_on_proc[get_from], MPI_INT,    win_numpar);
                    
            MPI_Get(&(remote_batches->ibeg[previous_remote_batches_counter]),   num_batches_on_proc[get_from], MPI_INT,
                    get_from, 0, num_batches_on_proc[get_from], MPI_INT,    win_ibeg);
                    
            MPI_Get(&(remote_batches->iend[previous_remote_batches_counter]),   num_batches_on_proc[get_from], MPI_INT,
                    get_from, 0, num_batches_on_proc[get_from], MPI_INT,    win_iend);

            MPI_Win_unlock(get_from, win_x_mid);
            MPI_Win_unlock(get_from, win_y_mid);
            MPI_Win_unlock(get_from, win_z_mid);
            MPI_Win_unlock(get_from, win_radius);
            MPI_Win_unlock(get_from, win_numpar);
            MPI_Win_unlock(get_from, win_ibeg);
            MPI_Win_unlock(get_from, win_iend);
            

            // Fill in LET tree array from Remote tree array.
            for (int i = 0; i < num_batches_on_proc[get_from]; ++i) {
                        
                // Set the beginning and ending particle indices for the associated nodes in the local sources list
                remote_batches->ibeg[previous_remote_batches_counter + i] += previous_remote_sources_counter;
                remote_batches->iend[previous_remote_batches_counter + i] += previous_remote_sources_counter;
            }
            
            previous_remote_sources_length[get_from] = previous_remote_sources_counter;
            
            previous_remote_batches_counter += num_batches_on_proc[get_from];
            previous_remote_sources_counter += num_sources_on_proc[get_from];

        } //end loop over num_procs


        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Win_free(&win_x_mid);
        MPI_Win_free(&win_y_mid);
        MPI_Win_free(&win_z_mid);
        MPI_Win_free(&win_radius);
        MPI_Win_free(&win_numpar);
        MPI_Win_free(&win_ibeg);
        MPI_Win_free(&win_iend);
        

        Particles_Alloc(&remote_sources, num_remote_sources);

        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

            int get_from = (num_procs+rank-proc_id) % num_procs;

            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_sources_x);
            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_sources_y);
            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_sources_z);
            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_sources_q);
            MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_sources_w);

            MPI_Get(&(remote_sources->x[previous_remote_sources_length[get_from]]),
                    num_sources_on_proc[get_from], MPI_DOUBLE,
                    get_from, 0, num_sources_on_proc[get_from], MPI_DOUBLE, win_sources_x);
            
            MPI_Get(&(remote_sources->y[previous_remote_sources_length[get_from]]),
                    num_sources_on_proc[get_from], MPI_DOUBLE,
                    get_from, 0, num_sources_on_proc[get_from], MPI_DOUBLE, win_sources_y);
            
            MPI_Get(&(remote_sources->z[previous_remote_sources_length[get_from]]),
                    num_sources_on_proc[get_from], MPI_DOUBLE,
                    get_from, 0, num_sources_on_proc[get_from], MPI_DOUBLE, win_sources_z);
            MPI_Get(&(remote_sources->q[previous_remote_sources_length[get_from]]),
                    num_sources_on_proc[get_from], MPI_DOUBLE,
                    get_from, 0, num_sources_on_proc[get_from], MPI_DOUBLE, win_sources_q);
            MPI_Get(&(remote_sources->w[previous_remote_sources_length[get_from]]),
                    num_sources_on_proc[get_from], MPI_DOUBLE,
                    get_from, 0, num_sources_on_proc[get_from], MPI_DOUBLE, win_sources_w);

            MPI_Win_unlock(get_from, win_sources_x);
            MPI_Win_unlock(get_from, win_sources_y);
            MPI_Win_unlock(get_from, win_sources_z);
            MPI_Win_unlock(get_from, win_sources_q);
            MPI_Win_unlock(get_from, win_sources_w);

        } // end loop over num_procs

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Win_free(&win_sources_x);
        MPI_Win_free(&win_sources_y);
        MPI_Win_free(&win_sources_z);
        MPI_Win_free(&win_sources_q);
        MPI_Win_free(&win_sources_w);

        time_tree[3] = MPI_Wtime() - time1;


        // Beginning local computation
        
        time1 = MPI_Wtime();

        struct InteractionLists *local_interaction_list = NULL;

        InteractionLists_Make(&local_interaction_list, tree_array, batches, run_params);

        time_tree[4] = MPI_Wtime() - time1; //time_constructlet


        time1 = MPI_Wtime();

        InteractionCompute_CP(tree_array, batches,
                        local_interaction_list,
                        sources->x, sources->y, sources->z, sources->q, sources->w,
                        targets->x, targets->y, targets->z, targets->q,
                        clusters->x, clusters->y, clusters->z, clusters->q, clusters->w,
                        potential_array, sources->num, targets->num, clusters->num,
                        run_params);
                        
        InteractionLists_Free(local_interaction_list);

        time_tree[5] = MPI_Wtime() - time1; //time_constructlet


        if (num_procs > 1) {
            time1 = MPI_Wtime();

            struct InteractionLists *let_interaction_list = NULL;
            
            InteractionLists_Make(&let_interaction_list, tree_array, remote_batches, run_params);

            time_tree[6] = MPI_Wtime() - time1; //time_makeglobintlist

            // After filling LET, call interaction_list_treecode
            time1 = MPI_Wtime(); // start timer for tree evaluation

            InteractionCompute_CP(tree_array, remote_batches,
                                    let_interaction_list,
                                    remote_sources->x, remote_sources->y, remote_sources->z, remote_sources->q, remote_sources->w,
                                    targets->x, targets->y, targets->z, targets->q,
                                    clusters->x, clusters->y, clusters->z, clusters->q, clusters->w,
                                    potential_array, remote_sources->num, targets->num, clusters->num,
                                    run_params);
            
            InteractionLists_Free(let_interaction_list);
        
            Particles_Free(remote_sources);
            Batches_Free(remote_batches);
            
            time_tree[7] = MPI_Wtime() - time1;
        }


        time1 = MPI_Wtime();

        InteractionCompute_Downpass(tree_array,
                                targets->x, targets->y, targets->z, targets->q,
                                clusters->x, clusters->y, clusters->z, clusters->q, clusters->w,
                                potential_array,
                                targets->num, clusters->num_charges, clusters->num_weights,
                                run_params);

        time_tree[8] = MPI_Wtime() - time1;

        time1 = MPI_Wtime();

        
        InteractionCompute_SubtractionPotentialCorrection(potential_array, targets->q, targets->num, run_params);
        Particles_ReorderTargetsAndPotential(targets, potential_array);

        time_tree[9] = MPI_Wtime() - time1;
        
        
        //-------------------
        //CLEANUP
        //-------------------
        

        time1 = MPI_Wtime();
        
        free_vector(targets->order); // free particle order arrays
        free_vector(sources->order); // free particle order arrays
        
        Tree_Free(troot);
        Tree_FreeArray(tree_array);
        
        //I'm still not sure about deallocating the clusters and batches here.
        //Clusters_Free_Win(clusters);
        //Batches_Free(batches);
        
        time_tree[10] = MPI_Wtime() - time1; //time_cleanup
        time_tree[11] = time_tree[0] + time_tree[1] + time_tree[3] + time_tree[4] + time_tree[6]; //total setup time
        time_tree[12] = time_tree[5] + time_tree[7] + time_tree[8]; // total compute time
    
        MPI_Barrier(MPI_COMM_WORLD);

        return;
        
        
        
        //--------------------------------------
        //--------------------------------------
        //PARTICLE CLUSTER
        //--------------------------------------
        //--------------------------------------


        
    } else if (run_params->compute_type == PARTICLE_CLUSTER) {
    
        struct tnode *troot = NULL;
        struct tnode_array *tree_array = NULL;
        double xyzminmax[6];
        int numnodes = 0, numleaves = 0;

        struct tnode_array *batches = NULL;
        double batch_lim[6];

        struct clusters *clusters = NULL;

        time1 = MPI_Wtime();

        Tree_Setup(sources, targets, run_params->interp_order, xyzminmax);
        Tree_PC_Create(&troot, sources, 1, sources->num,
                       run_params->max_per_source_leaf, xyzminmax, 0, &numnodes, &numleaves);
        Tree_SetIndex(troot, 0);
        Tree_AllocArray(&tree_array, numnodes);
        Tree_CreateArray(troot, tree_array);

        time_tree[0] = MPI_Wtime() - time1; //time_maketreearray
        

        time1 = MPI_Wtime();

        Batches_Alloc(&batches, batch_lim, targets, run_params->max_per_target_leaf);
        Batches_CreateTargetBatches(batches, targets, 1, targets->num, run_params->max_per_target_leaf, batch_lim);

        time_tree[1] = MPI_Wtime() - time1; //time_createbatch
        

        time1 = MPI_Wtime();

        Clusters_PC_Setup(&clusters, sources, run_params->interp_order, tree_array,
                          run_params->approximation, run_params->singularity);

        time_tree[2] = MPI_Wtime() - time1; //time_fillclusters


        //-------------------
        //BEGIN COMPUTE PHASE
        //-------------------

        MPI_Barrier(MPI_COMM_WORLD);
  
        time1 = MPI_Wtime();
        
        struct CommTypes *comm_types = NULL;
        struct CommWindows *comm_windows = NULL;
        struct tnode_array **let_tree_arrays = NULL;

        struct clusters *let_clusters = NULL;
        struct particles *let_sources = NULL;
        

        CommTypesAndTrees_Construct(&comm_types, &let_tree_arrays,
                                    tree_array, batches, run_params);


        Particles_Alloc(&let_sources, comm_types->let_sources_length);
        Clusters_Alloc(&let_clusters, comm_types->let_clusters_length,
                       run_params->approximation, run_params->singularity);
                                                

        CommWindows_Create(&comm_windows, clusters, sources);
        
        
        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

            int get_from = (num_procs + rank - proc_id) % num_procs;
            
            CommWindows_Lock(comm_windows, get_from);

            //This is a non-blocking call!
            CommWindows_GetData(let_clusters, let_sources, comm_types, comm_windows,
                                get_from, run_params);
                                
            CommWindows_Unlock(comm_windows, get_from);
        }
        
        CommWindows_Free(comm_windows);
        
        time_tree[3] = MPI_Wtime() - time1;


        //-------------------
        //LOCAL COMPUTE
        //-------------------
        
        time1 = MPI_Wtime();

        struct InteractionLists *local_interaction_list;
        
        InteractionLists_Make(&local_interaction_list, tree_array, batches, run_params);

        time_tree[4] = MPI_Wtime() - time1; //time_constructlet


        if (run_params->verbosity > 0) {
            for (int j = 0; j < batches->numnodes; j++){
                totalNumberApprox += local_interaction_list->num_approx[j];
                totalNumberDirect += local_interaction_list->num_direct[j];
            }
        }
        
        
        time1 = MPI_Wtime();

        InteractionCompute_PC(tree_array, batches,
                        local_interaction_list,
                        sources->x, sources->y, sources->z, sources->q, sources->w,
                        targets->x, targets->y, targets->z, targets->q,
                        clusters->x, clusters->y, clusters->z, clusters->q, clusters->w,
                        potential_array,
                        sources->num, targets->num, clusters->num,
                        run_params);
                        
        InteractionLists_Free(local_interaction_list);

        time_tree[5] = MPI_Wtime() - time1; //time_constructlet


        //-------------------
        //REMOTE COMPUTE
        //-------------------
        
        time_tree[6] = 0;
        time_tree[7] = 0;
            
        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {
            time1 = MPI_Wtime();
            int get_from = (num_procs+rank-proc_id) % num_procs;

            struct InteractionLists *let_interaction_list;
            InteractionLists_Make(&let_interaction_list, let_tree_arrays[get_from], batches, run_params);

            // Count number of interactions
            if (run_params->verbosity > 0) {
                for (int j = 0; j < batches->numnodes; j++) {
                    totalNumberApprox += let_interaction_list->num_approx[j];
                    totalNumberDirect += let_interaction_list->num_direct[j];
                }
            }
            time_tree[6] += MPI_Wtime() - time1; //time_makeglobintlist


            // After filling LET, call interaction_list_treecode
            time1 = MPI_Wtime(); // start timer for tree evaluation

            InteractionCompute_PC(let_tree_arrays[get_from], batches,
                                   let_interaction_list,
                                   let_sources->x, let_sources->y, let_sources->z, let_sources->q, let_sources->w,
                                   targets->x, targets->y, targets->z, targets->q,
                                   let_clusters->x, let_clusters->y, let_clusters->z, let_clusters->q, let_clusters->w,
                                   potential_array,
                                   let_sources->num, targets->num, let_clusters->num,
                                   run_params);
            
            InteractionLists_Free(let_interaction_list);
            
            time_tree[7] += MPI_Wtime() - time1;
            
        }
            
        
        //-------------------
        //CORRECT AND REORDER
        //-------------------

        time1 = MPI_Wtime();
        
        InteractionCompute_SubtractionPotentialCorrection(potential_array, targets->q, targets->num,
                                  run_params);

        Particles_ReorderTargetsAndPotential(targets, potential_array);

        time_tree[8] = 0.0;
        time_tree[9] = MPI_Wtime() - time1;
        
        
        //-------------------
        //CLEANUP
        //-------------------
        
        time1 = MPI_Wtime();
        
        free_vector(targets->order);
        free_vector(sources->order);
        
        Tree_Free(troot);
        Tree_FreeArray(tree_array);
        Clusters_Free_Win(clusters);
        Batches_Free(batches);
        
        //remote pieces
        Clusters_Free(let_clusters);
        Particles_Free(let_sources);
        CommTypesAndTrees_Free(comm_types, let_tree_arrays);
        
        
        time_tree[10] = MPI_Wtime() - time1; //time_cleanup
        time_tree[11] = time_tree[0] + time_tree[1] + time_tree[3] + time_tree[4] + time_tree[6]; //total setup time
        time_tree[12] = time_tree[5] + time_tree[7] + time_tree[8]; // total compute time
    
        MPI_Barrier(MPI_COMM_WORLD);

        return;



        //--------------------------------------
        //--------------------------------------
        //CLUSTER CLUSTER
        //--------------------------------------
        //--------------------------------------



    } else if (run_params->compute_type == CLUSTER_CLUSTER) {
    
        struct tnode *source_tree_root = NULL;
        struct tnode_array *source_tree_array = NULL;

        struct tnode *target_tree_root = NULL;
        struct tnode_array *target_tree_array = NULL;

        struct clusters *source_clusters = NULL;
        struct clusters *target_clusters = NULL;

        double source_xyzminmax[6], target_xyzminmax[6];
        int source_numnodes = 0, source_numleaves = 0;
        int target_numnodes = 0, target_numleaves = 0;

        time1 = MPI_Wtime();

        Tree_CC_Setup(sources, targets, run_params->interp_order,
                      source_xyzminmax, target_xyzminmax);

        Tree_PC_Create(&source_tree_root, sources, 1, sources->num,
                       run_params->max_per_source_leaf, source_xyzminmax, 0,
                       &source_numnodes, &source_numleaves);
        Tree_SetIndex(source_tree_root, 0);

        Tree_AllocArray(&source_tree_array, source_numnodes);
        Tree_CreateArray(source_tree_root, source_tree_array);

        time_tree[0] = MPI_Wtime() - time1;


        time1 = MPI_Wtime();

        Tree_CP_Create(&target_tree_root, targets, 1, targets->num,
                       run_params->max_per_target_leaf, target_xyzminmax, 0,
                       &target_numnodes, &target_numleaves);
        Tree_SetIndex(target_tree_root, 0);

        Tree_AllocArray(&target_tree_array, target_numnodes);
        Tree_CreateArray(target_tree_root, target_tree_array);

        time_tree[1] = MPI_Wtime() - time1; //time_maketreearray
         

        time1 = MPI_Wtime();

        Clusters_PC_Setup(&source_clusters, sources, run_params->interp_order, source_tree_array,
                          run_params->approximation, run_params->singularity);

        Clusters_CP_Setup(&target_clusters, run_params->interp_order, target_tree_array,
                          run_params->approximation, run_params->singularity);

        time_tree[2] = MPI_Wtime() - time1; //time_fillclusters
        
        
        //-------------------
        //COMPUTE PHASE
        //-------------------
        
        MPI_Barrier(MPI_COMM_WORLD);
    
        time1 = MPI_Wtime();
          
        struct CommTypes *comm_types = NULL;
        struct CommWindows *comm_windows = NULL;
        struct tnode_array **let_tree_arrays = NULL;

        struct clusters *let_clusters = NULL;
        struct particles *let_sources = NULL;

          
        CommTypesAndTrees_Construct(&comm_types, &let_tree_arrays,
                                    source_tree_array, target_tree_array, run_params);


        Particles_Alloc(&let_sources, comm_types->let_sources_length);
        Clusters_Alloc(&let_clusters, comm_types->let_clusters_length,
                       run_params->approximation, run_params->singularity);
                                                  

        CommWindows_Create(&comm_windows, source_clusters, sources);
        
        
        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

            int get_from = (num_procs + rank - proc_id) % num_procs;
            
            CommWindows_Lock(comm_windows, get_from);

            //This is a non-blocking call!
            CommWindows_GetData(let_clusters, let_sources, comm_types, comm_windows,
                                get_from, run_params);
                                
            CommWindows_Unlock(comm_windows, get_from);
        }
        
        CommWindows_Free(comm_windows);

        time_tree[3] = MPI_Wtime() - time1;

        
        //-------------------
        //LOCAL COMPUTE
        //-------------------
        
        time1 = MPI_Wtime();
        
        struct InteractionLists *local_interaction_list;
        InteractionLists_Make(&local_interaction_list, source_tree_array, target_tree_array, run_params);

        time_tree[4] = MPI_Wtime() - time1; //time_constructlet


        if (run_params->verbosity > 0) {
            for (int j = 0; j < target_tree_array->numnodes; j++){
                totalNumberApprox += local_interaction_list->num_approx[j];
                totalNumberDirect += local_interaction_list->num_direct[j];
            }
        }
        
        
        time1 = MPI_Wtime();

        InteractionCompute_CC(source_tree_array, target_tree_array,
                        local_interaction_list,
                        sources->x, sources->y, sources->z, sources->q, sources->w,
                        targets->x, targets->y, targets->z, targets->q,
                        source_clusters->x, source_clusters->y, source_clusters->z,
                        source_clusters->q, source_clusters->w,
                        target_clusters->x, target_clusters->y, target_clusters->z,
                        target_clusters->q, target_clusters->w,
                        potential_array,
                        sources->num, targets->num, source_clusters->num, target_clusters->num,
                        run_params);
                        
        InteractionLists_Free(local_interaction_list);

        time_tree[5] = MPI_Wtime() - time1; //time_constructlet
        

        //-------------------
        //REMOTE COMPUTE
        //-------------------
        
        time_tree[6] = 0;
        time_tree[7] = 0;
            
        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {
            time1 = MPI_Wtime();
            int get_from = (num_procs+rank-proc_id) % num_procs;

            struct InteractionLists *let_interaction_list;
            InteractionLists_Make(&let_interaction_list, let_tree_arrays[get_from], target_tree_array,  run_params);

            time_tree[6] += MPI_Wtime() - time1; //time_makeglobintlist
             
             
            // Count number of interactions
            if (run_params->verbosity > 0) {
                for (int j = 0; j < target_tree_array->numnodes; j++){
                    totalNumberApprox += let_interaction_list->num_approx[j];
                    totalNumberDirect += let_interaction_list->num_direct[j];
                }
            }
           

            // After filling LET, call interaction_list_treecode
            time1 = MPI_Wtime(); // start timer for tree evaluation

            InteractionCompute_CC(let_tree_arrays[get_from], target_tree_array,
                                   let_interaction_list,
                                   let_sources->x, let_sources->y, let_sources->z,
                                   let_sources->q, let_sources->w,
                                   targets->x, targets->y, targets->z, targets->q,
                                   let_clusters->x, let_clusters->y, let_clusters->z,
                                   let_clusters->q, let_clusters->w,
                                   target_clusters->x, target_clusters->y, target_clusters->z,
                                   target_clusters->q, target_clusters->w,
                                   potential_array,
                                   let_sources->num, targets->num, let_clusters->num, target_clusters->num,
                                   run_params);

            InteractionLists_Free(let_interaction_list);
            
            time_tree[7] += MPI_Wtime() - time1;
            
        }
            
        
        //-------------------
        //DOWNPASS
        //-------------------
        
        time1 = MPI_Wtime();
        
        InteractionCompute_Downpass(target_tree_array,
                                targets->x, targets->y, targets->z, targets->q,
                                target_clusters->x, target_clusters->y, target_clusters->z,
                                target_clusters->q, target_clusters->w,
                                potential_array, 
                                targets->num, target_clusters->num_charges, target_clusters->num_weights,
                                run_params);

        time_tree[8] = MPI_Wtime() - time1;
            
        
        //-------------------
        //CORRECT AND REORDER
        //-------------------

        time1 = MPI_Wtime();
        
        InteractionCompute_SubtractionPotentialCorrection(potential_array, targets->q, targets->num,
                                  run_params);

        Particles_ReorderTargetsAndPotential(targets, potential_array);

        time_tree[9] = MPI_Wtime() - time1;
               
               
        //-------------------
        //CLEANUP
        //-------------------
        
        time1 = MPI_Wtime();
        
        free_vector(targets->order);
        free_vector(sources->order);

        Tree_Free(source_tree_root);
        Tree_Free(target_tree_root);
        Tree_FreeArray(source_tree_array);
        Tree_FreeArray(target_tree_array);
        
        Clusters_Free_Win(source_clusters);
        Clusters_Free(target_clusters);
        
        //remote pieces
        Clusters_Free(let_clusters);
        Particles_Free(let_sources);
        CommTypesAndTrees_Free(comm_types, let_tree_arrays);
        
        time_tree[10] = MPI_Wtime() - time1; //time_cleanup
        time_tree[11] = time_tree[0] + time_tree[1] + time_tree[3] + time_tree[4] + time_tree[6]; //total setup time
        time_tree[12] = time_tree[5] + time_tree[7] + time_tree[8]; // total compute time
        
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
//        totalNumberInteractions=totalNumberDirect+totalNumberApprox;
//        printf("Interaction information: \n");
//        printf("rank %d: number of direct batch-cluster interactions: %d\n", rank, totalNumberApprox);
//        printf("rank %d: number of approx batch-cluster interactions: %d\n", rank, totalNumberDirect);
//        printf("rank %d:  total number of batch-cluster interactions: %d\n\n", rank, totalNumberInteractions);
//        MPI_Barrier(MPI_COMM_WORLD);
//    }
//
//    if (run_params -> verbosity > 0) {
//        MPI_Reduce(&totalNumberInteractions,&cumulativeNumberInteractions, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
//        MPI_Reduce(&totalNumberInteractions,&maxNumberInteractions, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
//        MPI_Reduce(&totalNumberInteractions,&minNumberInteractions, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
//        if (rank==0){
//            printf("Cumulative number of interactions across all ranks: %d\n", cumulativeNumberInteractions);
//            printf("   Maximum number of interactions across all ranks: %d\n", maxNumberInteractions);
//            printf("   Minimum number of interactions across all ranks: %d\n", minNumberInteractions);
//            printf("                                             Ratio: %f\n\n", (double)maxNumberInteractions/(double)minNumberInteractions );
//        }
//        MPI_Barrier(MPI_COMM_WORLD);
//    }

} /* END function treecode */
