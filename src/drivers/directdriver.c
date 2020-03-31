#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../utilities/timers.h"
#include "../utilities/enums.h"

#include "../particles/struct_particles.h"
#include "../particles/particles.h"

#include "../run_params/struct_run_params.h"

#include "../interaction_compute/interaction_compute.h"

#include "directdriver.h"


void directdriver(struct Particles *sources, struct Particles *targets, struct RunParams *run_params,
                  double *potential, double *time_direct)
{
    int rank, num_procs, ierr;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int num_sources = sources->num;
    int num_targets = targets->num;
    int num_sources_on_proc[num_procs];

    double *source_x = sources->x;
    double *source_y = sources->y;
    double *source_z = sources->z;
    double *source_q = sources->q;
    double *source_w = sources->w;

    double *target_x = targets->x;
    double *target_y = targets->y;
    double *target_z = targets->z;
    double *target_q = targets->q;

    double time1;
    time_direct[0] = 0.0;
    time_direct[1] = 0.0;
    time_direct[2] = 0.0;
    time_direct[3] = 0.0;


    START_TIMER(&time1);
    MPI_Allgather(&num_sources, 1, MPI_INT, num_sources_on_proc, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Win win_sources_x, win_sources_y, win_sources_z, win_sources_q, win_sources_w;
    MPI_Win_create(source_x, num_sources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_x);
    MPI_Win_create(source_y, num_sources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_y);
    MPI_Win_create(source_z, num_sources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_z);
    MPI_Win_create(source_q, num_sources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_q);
    MPI_Win_create(source_w, num_sources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_w);
    STOP_TIMER(&time1);
    time_direct[0] += time1;

    for (int proc_id = 1; proc_id < num_procs; ++proc_id) {

        START_TIMER(&time1);
        int get_from = (num_procs + rank - proc_id) % num_procs;

        struct Particles *remote_sources = NULL;
        Particles_Alloc(&remote_sources, num_sources_on_proc[get_from]);
        
        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_sources_x);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_sources_y);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_sources_z);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_sources_q);
        MPI_Win_lock(MPI_LOCK_SHARED, get_from, 0, win_sources_w);

        MPI_Get(remote_sources->x, num_sources_on_proc[get_from], MPI_DOUBLE,
                       get_from, 0, num_sources_on_proc[get_from], MPI_DOUBLE, win_sources_x);
        MPI_Get(remote_sources->y, num_sources_on_proc[get_from], MPI_DOUBLE,
                       get_from, 0, num_sources_on_proc[get_from], MPI_DOUBLE, win_sources_y);
        MPI_Get(remote_sources->z, num_sources_on_proc[get_from], MPI_DOUBLE,
                       get_from, 0, num_sources_on_proc[get_from], MPI_DOUBLE, win_sources_z);
        MPI_Get(remote_sources->q, num_sources_on_proc[get_from], MPI_DOUBLE,
                       get_from, 0, num_sources_on_proc[get_from], MPI_DOUBLE, win_sources_q);
        MPI_Get(remote_sources->w, num_sources_on_proc[get_from], MPI_DOUBLE,
                       get_from, 0, num_sources_on_proc[get_from], MPI_DOUBLE, win_sources_w);

        MPI_Win_unlock(get_from, win_sources_x);
        MPI_Win_unlock(get_from, win_sources_y);
        MPI_Win_unlock(get_from, win_sources_z);
        MPI_Win_unlock(get_from, win_sources_q);
        MPI_Win_unlock(get_from, win_sources_w);
        
        MPI_Barrier(MPI_COMM_WORLD);

        STOP_TIMER(&time1);
        time_direct[0] += time1;


        START_TIMER(&time1);
        InteractionCompute_Direct(potential, remote_sources, targets, run_params);

        Particles_Free(&remote_sources);
        STOP_TIMER(&time1);
        time_direct[1] += time1;
    }


    START_TIMER(&time_direct[2]);
    InteractionCompute_Direct(potential, sources, targets, run_params);
    STOP_TIMER(&time_direct[2]);


    START_TIMER(&time_direct[3]);
    InteractionCompute_SubtractionPotentialCorrection(potential, targets, run_params);
    STOP_TIMER(&time_direct[3]);

    return;

} /* END of function directdriver */
