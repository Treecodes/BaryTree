#include <stdio.h>
#include <mpi.h>

#include "../utilities/array.h"
#include "../tree/struct_tree.h"
#include "../tree/tree.h"
#include "../interaction_lists/interaction_lists.h"
#include "../run_params/struct_run_params.h"

#include "struct_comm_types.h"


void CommTypesAndTrees_Construct(struct CommTypes **comm_types_addr, struct Tree ***let_trees_addr,
                                 struct Tree *tree, struct Tree *batches,
                                 struct RunParams *run_params)
{
    MPI_Barrier(MPI_COMM_WORLD);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
        
    int *num_nodes_on_proc;
    make_vector(num_nodes_on_proc, num_procs);
    MPI_Allgather(&(tree->numnodes), 1, MPI_INT, num_nodes_on_proc, 1, MPI_INT, MPI_COMM_WORLD);


    MPI_Win win_x_mid, win_y_mid, win_z_mid, win_radius, win_numpar, win_ibeg, win_iend;
    MPI_Win win_children, win_num_children;


    *comm_types_addr = malloc(sizeof (struct CommTypes));  
    struct CommTypes *comm_types = *comm_types_addr;
    
    comm_types->let_clusters_length = 0;
    comm_types->let_clusters_num = 0;
    comm_types->let_sources_length = 0;

    make_vector(comm_types->num_remote_approx_array, num_procs);
    make_vector(comm_types->previous_let_clusters_length_array, num_procs);

    make_vector(comm_types->MPI_approx_type, num_procs);
    make_vector(comm_types->MPI_approx_charges_type, num_procs);
    make_vector(comm_types->MPI_approx_weights_type, num_procs);

    make_vector(comm_types->new_sources_length_array, num_procs);
    make_vector(comm_types->previous_let_sources_length_array, num_procs);

    make_vector(comm_types->MPI_direct_type, num_procs);

    
    *let_trees_addr = malloc(num_procs * sizeof (struct Tree *));
    struct Tree **let_trees = *let_trees_addr;


    MPI_Win_create(tree->x_mid,         tree->numnodes*sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_x_mid);

    MPI_Win_create(tree->y_mid,         tree->numnodes*sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_y_mid);
    
    MPI_Win_create(tree->z_mid,         tree->numnodes*sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_z_mid);
    
    MPI_Win_create(tree->radius,        tree->numnodes*sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_radius);
    
    MPI_Win_create(tree->numpar,        tree->numnodes*sizeof(int),    sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_numpar);
    
    MPI_Win_create(tree->ibeg,          tree->numnodes*sizeof(int),    sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_ibeg);
    
    MPI_Win_create(tree->iend,          tree->numnodes*sizeof(int),    sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_iend);
    
    MPI_Win_create(tree->num_children,  tree->numnodes*sizeof(int),    sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_num_children);
    
    MPI_Win_create(tree->children,    8*tree->numnodes*sizeof(int),    sizeof(int),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &win_children);


    for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

        int get_from = (num_procs + rank - proc_id) % num_procs;
        
        let_trees[get_from] = NULL;
        Tree_Alloc(&(let_trees[get_from]), num_nodes_on_proc[get_from]);
        struct Tree *remote_tree = let_trees[get_from];

        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_x_mid);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_y_mid);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_z_mid);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_radius);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_numpar);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_ibeg);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_iend);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_children);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_num_children);
        
        MPI_Get(remote_tree->x_mid,        num_nodes_on_proc[get_from], MPI_DOUBLE,
                get_from, 0, num_nodes_on_proc[get_from],   MPI_DOUBLE, win_x_mid);

        MPI_Get(remote_tree->y_mid,        num_nodes_on_proc[get_from], MPI_DOUBLE,
                get_from, 0, num_nodes_on_proc[get_from],   MPI_DOUBLE, win_y_mid);

        MPI_Get(remote_tree->z_mid,        num_nodes_on_proc[get_from], MPI_DOUBLE,
                get_from, 0, num_nodes_on_proc[get_from],   MPI_DOUBLE, win_z_mid);

        MPI_Get(remote_tree->radius,       num_nodes_on_proc[get_from], MPI_DOUBLE,
                get_from, 0, num_nodes_on_proc[get_from],   MPI_DOUBLE, win_radius);

        MPI_Get(remote_tree->numpar,       num_nodes_on_proc[get_from], MPI_INT,
                get_from, 0, num_nodes_on_proc[get_from],   MPI_INT,    win_numpar);

        MPI_Get(remote_tree->ibeg,         num_nodes_on_proc[get_from], MPI_INT,
                get_from, 0, num_nodes_on_proc[get_from],   MPI_INT,    win_ibeg);

        MPI_Get(remote_tree->iend,         num_nodes_on_proc[get_from], MPI_INT,
                get_from, 0, num_nodes_on_proc[get_from],   MPI_INT,    win_iend);

        MPI_Get(remote_tree->children,   8*num_nodes_on_proc[get_from], MPI_INT,
                get_from, 0, 8*num_nodes_on_proc[get_from], MPI_INT,    win_children);

        MPI_Get(remote_tree->num_children, num_nodes_on_proc[get_from], MPI_INT,
                get_from, 0, num_nodes_on_proc[get_from],   MPI_INT,    win_num_children);

        MPI_Win_unlock(get_from, win_x_mid);
        MPI_Win_unlock(get_from, win_y_mid);
        MPI_Win_unlock(get_from, win_z_mid);
        MPI_Win_unlock(get_from, win_radius);
        MPI_Win_unlock(get_from, win_numpar);
        MPI_Win_unlock(get_from, win_ibeg);
        MPI_Win_unlock(get_from, win_iend);
        MPI_Win_unlock(get_from, win_children);
        MPI_Win_unlock(get_from, win_num_children);


        int *approx_list_packed, *approx_list_unpacked, *direct_list, *direct_ibeg_list, *direct_length_list;
        make_vector(approx_list_packed, num_nodes_on_proc[get_from]);
        make_vector(approx_list_unpacked, num_nodes_on_proc[get_from]);
        make_vector(direct_list, num_nodes_on_proc[get_from]);
        make_vector(direct_ibeg_list, num_nodes_on_proc[get_from]);
        make_vector(direct_length_list, num_nodes_on_proc[get_from]);


        InteractionLists_MakeRemote(remote_tree, batches, approx_list_unpacked, approx_list_packed,
                                    direct_list, run_params);


        int num_remote_approx = 0;
        int previous_let_clusters_length = comm_types->let_clusters_length;

        int num_remote_direct = 0;
        int previous_let_sources_length = comm_types->let_sources_length;
        
        
        int append_counter = 0;
        for (int i = 0; i < num_nodes_on_proc[get_from]; ++i) {

            if (approx_list_unpacked[i] != -1) {
            
                remote_tree->cluster_ind[i] = comm_types->let_clusters_num;
                comm_types->let_clusters_length += run_params->interp_pts_per_cluster;
                comm_types->let_clusters_num++;
                num_remote_approx++;
            }
                     
            if (direct_list[i] != -1) {
                 
                // Determine displacements and lengths for getting particles from remote sources list
                direct_ibeg_list[num_remote_direct] = remote_tree->ibeg[i] - 1; // zero index based
                direct_length_list[num_remote_direct] = remote_tree->numpar[i];
                num_remote_direct++;
        
                // Set beginning and ending particle indices for associated nodes in let sources list
                remote_tree->ibeg[i] = comm_types->let_sources_length + 1; //one index based, for some reason
                remote_tree->iend[i] = comm_types->let_sources_length + remote_tree->numpar[i];
                comm_types->let_sources_length += remote_tree->numpar[i];
            }
                 
            append_counter++;
        }
        
        
        comm_types->num_remote_approx_array[get_from] = num_remote_approx;
        comm_types->new_sources_length_array[get_from] = comm_types->let_sources_length - previous_let_sources_length;
        comm_types->previous_let_clusters_length_array[get_from] = previous_let_clusters_length;
        comm_types->previous_let_sources_length_array[get_from] = previous_let_sources_length;
        
        int *approx_list_displacements, *approx_charges_list_displacements, *approx_weights_list_displacements;
        make_vector(approx_list_displacements, num_nodes_on_proc[get_from]);
        make_vector(approx_charges_list_displacements, num_nodes_on_proc[get_from]);
        make_vector(approx_weights_list_displacements, num_nodes_on_proc[get_from]);

        // Use masks to get remote data
        for (int ii = 0; ii < num_remote_approx; ++ii) {
            approx_list_displacements[ii] = approx_list_packed[ii] * run_params->interp_pts_per_cluster;
            approx_charges_list_displacements[ii] = approx_list_packed[ii] * run_params->interp_charges_per_cluster;
            approx_weights_list_displacements[ii] = approx_list_packed[ii] * run_params->interp_weights_per_cluster;
        }
        
        MPI_Type_create_indexed_block(num_remote_approx, run_params->interp_pts_per_cluster,
                                      approx_list_displacements,
                                      MPI_DOUBLE, &(comm_types->MPI_approx_type[get_from]));
        MPI_Type_commit(&(comm_types->MPI_approx_type[get_from]));

        MPI_Type_create_indexed_block(num_remote_approx, run_params->interp_charges_per_cluster,
                                      approx_charges_list_displacements,
                                      MPI_DOUBLE, &(comm_types->MPI_approx_charges_type[get_from]));
        MPI_Type_commit(&(comm_types->MPI_approx_charges_type[get_from]));

        MPI_Type_create_indexed_block(num_remote_approx, run_params->interp_weights_per_cluster,
                                      approx_weights_list_displacements,
                                      MPI_DOUBLE, &(comm_types->MPI_approx_weights_type[get_from]));
        MPI_Type_commit(&(comm_types->MPI_approx_weights_type[get_from]));

        MPI_Type_indexed(num_remote_direct, direct_length_list, direct_ibeg_list,
                                      MPI_DOUBLE, &(comm_types->MPI_direct_type[get_from]));
        MPI_Type_commit(&(comm_types->MPI_direct_type[get_from]));

        free_vector(approx_list_packed);
        free_vector(approx_list_unpacked);
        free_vector(approx_list_displacements);
        free_vector(approx_charges_list_displacements);
        free_vector(approx_weights_list_displacements);
        free_vector(direct_list);
        free_vector(direct_ibeg_list);
        free_vector(direct_length_list);
    } //end loop over numProcs
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Win_free(&win_x_mid);
    MPI_Win_free(&win_y_mid);
    MPI_Win_free(&win_z_mid);
    MPI_Win_free(&win_radius);
    MPI_Win_free(&win_numpar);
    MPI_Win_free(&win_ibeg);
    MPI_Win_free(&win_iend);
    MPI_Win_free(&win_children);
    MPI_Win_free(&win_num_children);

    return;
}



void CommTypesAndTrees_Free(struct CommTypes **comm_types_addr, struct Tree ***let_trees_addr)
{
    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    
    struct CommTypes *comm_types = *comm_types_addr;
    struct Tree **let_trees = *let_trees_addr;

    if (comm_types != NULL) {
        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

            int get_from = (num_procs + rank - proc_id) % num_procs;

            MPI_Type_free(&(comm_types->MPI_approx_type[get_from]));
            MPI_Type_free(&(comm_types->MPI_approx_charges_type[get_from]));
            MPI_Type_free(&(comm_types->MPI_approx_weights_type[get_from]));
            MPI_Type_free(&(comm_types->MPI_direct_type[get_from]));
        }
        
        free_vector(comm_types->num_remote_approx_array);
        free_vector(comm_types->previous_let_clusters_length_array);
        free_vector(comm_types->MPI_approx_type);
        free_vector(comm_types->MPI_approx_charges_type);
        free_vector(comm_types->MPI_approx_weights_type);
        free_vector(comm_types->new_sources_length_array);
        free_vector(comm_types->previous_let_sources_length_array);
        free_vector(comm_types->MPI_direct_type);
        
        free(comm_types);
        comm_types = NULL;
    }
    
    if (let_trees != NULL) {
        for (int proc_id = 1; proc_id < num_procs; ++proc_id) {
            Tree_Free(&let_trees[(num_procs + rank - proc_id) % num_procs]);
        }
        
        free_vector(let_trees);
        let_trees = NULL;
    }

    return;
}
