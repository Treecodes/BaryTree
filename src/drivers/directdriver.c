#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "../utilities/array.h"
#include "../utilities/tools.h"
#include "../utilities/enums.h"

#include "../particles/struct_particles.h"
#include "../particles/particles.h"

#include "../run_params/struct_run_params.h"

#include "../interaction_compute/interaction_compute.h"

#include "directdriver.h"


void directdriver(struct particles *sources, struct particles *targets, struct RunParams *run_params,
                  double *pointwisePotential, double *time_direct)
{
    int rank, numProcs, ierr;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int numSources = sources->num;
    int numTargets = targets->num;
    int numSourcesOnProc[numProcs];

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


    time1 = MPI_Wtime();

    MPI_Allgather(&numSources, 1, MPI_INT, numSourcesOnProc, 1, MPI_INT, MPI_COMM_WORLD);
    MPI_Win win_sources_x, win_sources_y, win_sources_z, win_sources_q, win_sources_w;
    MPI_Win_create(source_x, numSources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_x);
    MPI_Win_create(source_y, numSources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_y);
    MPI_Win_create(source_z, numSources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_z);
    MPI_Win_create(source_q, numSources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_q);
    MPI_Win_create(source_w, numSources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_w);

    time_direct[0] += MPI_Wtime() - time1;


    for (int procID = 1; procID < numProcs; ++procID) {

        time1 = MPI_Wtime();

        int getFrom = (numProcs+rank-procID) % numProcs;

        struct particles *remote_sources = malloc(sizeof(struct particles));
        Particles_AllocSources(remote_sources, numSourcesOnProc[getFrom]);

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_sources_x);
        MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_sources_y);
        MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_sources_z);
        MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_sources_q);
        MPI_Win_lock(MPI_LOCK_SHARED, getFrom, 0, win_sources_w);

        MPI_Get(remote_sources->x, numSourcesOnProc[getFrom], MPI_DOUBLE,
                       getFrom, 0, numSourcesOnProc[getFrom], MPI_DOUBLE, win_sources_x);
        MPI_Get(remote_sources->y, numSourcesOnProc[getFrom], MPI_DOUBLE,
                       getFrom, 0, numSourcesOnProc[getFrom], MPI_DOUBLE, win_sources_y);
        MPI_Get(remote_sources->z, numSourcesOnProc[getFrom], MPI_DOUBLE,
                       getFrom, 0, numSourcesOnProc[getFrom], MPI_DOUBLE, win_sources_z);
        MPI_Get(remote_sources->q, numSourcesOnProc[getFrom], MPI_DOUBLE,
                       getFrom, 0, numSourcesOnProc[getFrom], MPI_DOUBLE, win_sources_q);
        MPI_Get(remote_sources->w, numSourcesOnProc[getFrom], MPI_DOUBLE,
                       getFrom, 0, numSourcesOnProc[getFrom], MPI_DOUBLE, win_sources_w);

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Win_unlock(getFrom, win_sources_x);
        MPI_Win_unlock(getFrom, win_sources_y);
        MPI_Win_unlock(getFrom, win_sources_z);
        MPI_Win_unlock(getFrom, win_sources_q);
        MPI_Win_unlock(getFrom, win_sources_w);

        MPI_Barrier(MPI_COMM_WORLD);

        time_direct[0] += MPI_Wtime() - time1;

        
        time1 = MPI_Wtime();
        //compute remote

        fprintf(stderr, "I shouldn't be in here if I'm running in serial.\n");

        InteractionCompute_Direct(remote_sources->x, remote_sources->y, remote_sources->z,
                                   remote_sources->q, remote_sources->w,
                                   target_x, target_y, target_z, target_q,
                                   pointwisePotential, numSources, numTargets,
                                   run_params);

        Particles_FreeSources(remote_sources);

        time_direct[1] += MPI_Wtime() - time1;
    }


    time1 = MPI_Wtime();
    //compute local 
    InteractionCompute_Direct(source_x, source_y, source_z, source_q, source_w,
                               target_x, target_y, target_z, target_q,
                               pointwisePotential, numSources, numTargets,
                               run_params);

    time_direct[2] = MPI_Wtime() - time1;


    time1 = MPI_Wtime();
    //add correction
    InteractionCompute_SubtractionPotentialCorrection(pointwisePotential, target_q, numTargets, run_params);

    time_direct[3] = MPI_Wtime() - time1;

    return;

} /* END of function directdriver */
