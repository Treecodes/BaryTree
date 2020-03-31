#include <stdio.h>
#include <mpi.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../tree/struct_tree.h"
#include "../tree/batches.h"
#include "../particles/struct_particles.h"
#include "../particles/particles.h"
#include "../run_params/struct_run_params.h"


void Comm_CP_ConstructAndGetData(struct Tree **remote_batches_addr, struct Particles **remote_sources_addr,
                                 const struct Tree *tree, const struct Tree *batches,
                                 const struct Particles *sources, const struct RunParams *run_params)
{
    MPI_Barrier(MPI_COMM_WORLD);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
        
    int *num_batches_on_proc;
    make_vector(num_batches_on_proc, num_procs);
    MPI_Allgather(&(batches->numnodes), 1, MPI_INT, num_batches_on_proc, 1, MPI_INT, MPI_COMM_WORLD);

    int *num_sources_on_proc;
    make_vector(num_sources_on_proc, num_procs);
    MPI_Allgather(&(sources->num), 1, MPI_INT, num_sources_on_proc, 1, MPI_INT, MPI_COMM_WORLD);

    int *previous_remote_sources_length;
    make_vector(previous_remote_sources_length, num_procs);
    
    
    int num_remote_batches = sum_int(num_batches_on_proc, num_procs) - num_batches_on_proc[rank];
    int num_remote_sources = sum_int(num_sources_on_proc, num_procs) - num_sources_on_proc[rank];
    
    
    Batches_Alloc(remote_batches_addr, num_remote_batches);
    struct Tree *remote_batches = *remote_batches_addr;
    
    Particles_Alloc(remote_sources_addr, num_remote_sources);
    struct Particles *remote_sources = *remote_sources_addr;
    
    
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
        

        for (int i = 0; i < num_batches_on_proc[get_from]; ++i) {
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

    
    
    return;
}
