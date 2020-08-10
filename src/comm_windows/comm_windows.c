#include <stdlib.h>
#include <mpi.h>

#include "../clusters/struct_clusters.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"
#include "../comm_types/struct_comm_types.h"

#include "struct_comm_windows.h"

void CommWindows_Create(struct CommWindows **comm_windows_addr,
                        struct Clusters *clusters, struct Particles *sources, struct RunParams *run_params)
{
    *comm_windows_addr = malloc(sizeof(struct CommWindows));
    struct CommWindows *comm_windows = *comm_windows_addr;

    MPI_Win_create(clusters->x, clusters->num         * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &(comm_windows->win_clusters_x));

    MPI_Win_create(clusters->y, clusters->num         * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &(comm_windows->win_clusters_y));
    
    MPI_Win_create(clusters->z, clusters->num         * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &(comm_windows->win_clusters_z));
    
    MPI_Win_create(clusters->q, clusters->num_charges * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &(comm_windows->win_clusters_q));
    
    MPI_Win_create(sources->x,  sources->num          * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &(comm_windows->win_sources_x));
    
    MPI_Win_create(sources->y,  sources->num          * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &(comm_windows->win_sources_y));
    
    MPI_Win_create(sources->z,  sources->num          * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &(comm_windows->win_sources_z));
    
    MPI_Win_create(sources->q,  sources->num          * sizeof(double), sizeof(double),
                   MPI_INFO_NULL, MPI_COMM_WORLD, &(comm_windows->win_sources_q));
    
    if (run_params->singularity == SUBTRACTION) {
        MPI_Win_create(clusters->w, clusters->num_weights * sizeof(double), sizeof(double),
                       MPI_INFO_NULL, MPI_COMM_WORLD, &(comm_windows->win_clusters_w));

        MPI_Win_create(sources->w,  sources->num          * sizeof(double), sizeof(double),
                       MPI_INFO_NULL, MPI_COMM_WORLD, &(comm_windows->win_sources_w));
    }

    return;
}



void CommWindows_Free(struct CommWindows **comm_windows_addr, struct RunParams *run_params)
{
    MPI_Barrier(MPI_COMM_WORLD);
    struct CommWindows *comm_windows = *comm_windows_addr;

    MPI_Win_free(&(comm_windows->win_clusters_x));
    MPI_Win_free(&(comm_windows->win_clusters_y));
    MPI_Win_free(&(comm_windows->win_clusters_z));
    MPI_Win_free(&(comm_windows->win_clusters_q));

    MPI_Win_free(&(comm_windows->win_sources_x));
    MPI_Win_free(&(comm_windows->win_sources_y));
    MPI_Win_free(&(comm_windows->win_sources_z));
    MPI_Win_free(&(comm_windows->win_sources_q));

    if (run_params->singularity == SUBTRACTION) {
        MPI_Win_free(&(comm_windows->win_clusters_w));
        MPI_Win_free(&(comm_windows->win_sources_w));
    }

    free(comm_windows);
    comm_windows = NULL;

    return;
}



void CommWindows_Lock(struct CommWindows *comm_windows, int get_from, struct RunParams *run_params)
{
    MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, comm_windows->win_clusters_x);
    MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, comm_windows->win_clusters_y);
    MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, comm_windows->win_clusters_z);
    MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, comm_windows->win_clusters_q);

    MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, comm_windows->win_sources_x);
    MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, comm_windows->win_sources_y);
    MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, comm_windows->win_sources_z);
    MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, comm_windows->win_sources_q);

    if (run_params->singularity == SUBTRACTION) {
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, comm_windows->win_clusters_w);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, comm_windows->win_sources_w);
    }

    return;
}



void CommWindows_Unlock(struct CommWindows *comm_windows, int get_from, struct RunParams *run_params)
{
    MPI_Win_unlock(get_from, comm_windows->win_clusters_x);
    MPI_Win_unlock(get_from, comm_windows->win_clusters_y);
    MPI_Win_unlock(get_from, comm_windows->win_clusters_z);
    MPI_Win_unlock(get_from, comm_windows->win_clusters_q);

    MPI_Win_unlock(get_from, comm_windows->win_sources_x);
    MPI_Win_unlock(get_from, comm_windows->win_sources_y);
    MPI_Win_unlock(get_from, comm_windows->win_sources_z);
    MPI_Win_unlock(get_from, comm_windows->win_sources_q);

    if (run_params->singularity == SUBTRACTION) {
        MPI_Win_unlock(get_from, comm_windows->win_clusters_w);
        MPI_Win_unlock(get_from, comm_windows->win_sources_w);
    }

    return;
}



void CommWindows_GetData(struct Clusters *let_clusters, struct Particles *let_sources,
                         struct CommTypes *comm_types, struct CommWindows *comm_windows,
                         int get_from, struct RunParams *run_params)
{
    int interp_pts_per_cluster     = run_params->interp_pts_per_cluster;
    int interp_weights_per_cluster = run_params->interp_weights_per_cluster;
    int interp_charges_per_cluster = run_params->interp_charges_per_cluster;

    int weights_per_point = interp_weights_per_cluster / interp_pts_per_cluster;
    int charges_per_point = interp_charges_per_cluster / interp_pts_per_cluster;


    MPI_Get(&(let_clusters->x[comm_types->previous_let_clusters_length_array[get_from]]),
            comm_types->num_remote_approx_array[get_from] * interp_pts_per_cluster,     MPI_DOUBLE,
            get_from, 0, 1, comm_types->MPI_approx_type[get_from],         comm_windows->win_clusters_x);

    MPI_Get(&(let_clusters->y[comm_types->previous_let_clusters_length_array[get_from]]),
            comm_types->num_remote_approx_array[get_from] * interp_pts_per_cluster,     MPI_DOUBLE,
            get_from, 0, 1, comm_types->MPI_approx_type[get_from],         comm_windows->win_clusters_y);

    MPI_Get(&(let_clusters->z[comm_types->previous_let_clusters_length_array[get_from]]),
            comm_types->num_remote_approx_array[get_from] * interp_pts_per_cluster,     MPI_DOUBLE,
            get_from, 0, 1, comm_types->MPI_approx_type[get_from],         comm_windows->win_clusters_z);

    MPI_Get(&(let_clusters->q[comm_types->previous_let_clusters_length_array[get_from] * charges_per_point]),
            comm_types->num_remote_approx_array[get_from] * interp_charges_per_cluster, MPI_DOUBLE,
            get_from, 0, 1, comm_types->MPI_approx_charges_type[get_from], comm_windows->win_clusters_q);


    MPI_Get(&(let_sources->x[comm_types->previous_let_sources_length_array[get_from]]),
            comm_types->new_sources_length_array[get_from], MPI_DOUBLE,
            get_from, 0, 1, comm_types->MPI_direct_type[get_from], comm_windows->win_sources_x);

    MPI_Get(&(let_sources->y[comm_types->previous_let_sources_length_array[get_from]]),
            comm_types->new_sources_length_array[get_from], MPI_DOUBLE,
            get_from, 0, 1, comm_types->MPI_direct_type[get_from], comm_windows->win_sources_y);

    MPI_Get(&(let_sources->z[comm_types->previous_let_sources_length_array[get_from]]),
            comm_types->new_sources_length_array[get_from], MPI_DOUBLE,
            get_from, 0, 1, comm_types->MPI_direct_type[get_from], comm_windows->win_sources_z);

    MPI_Get(&(let_sources->q[comm_types->previous_let_sources_length_array[get_from]]),
            comm_types->new_sources_length_array[get_from], MPI_DOUBLE,
            get_from, 0, 1, comm_types->MPI_direct_type[get_from], comm_windows->win_sources_q);


    if (run_params->singularity == SUBTRACTION) {
        MPI_Get(&(let_clusters->w[comm_types->previous_let_clusters_length_array[get_from] * weights_per_point]),
                comm_types->num_remote_approx_array[get_from] * interp_weights_per_cluster, MPI_DOUBLE,
                get_from, 0, 1, comm_types->MPI_approx_weights_type[get_from], comm_windows->win_clusters_w);

        MPI_Get(&(let_sources->w[comm_types->previous_let_sources_length_array[get_from]]),
                comm_types->new_sources_length_array[get_from], MPI_DOUBLE,
                get_from, 0, 1, comm_types->MPI_direct_type[get_from], comm_windows->win_sources_w);
    }

    return;
}
