/*
 *Procedures for Particle-Cluster Treecode
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <mpi.h>

#include "array.h"
#include "tools.h"
#include "struct_particles.h"

#include "interaction_compute.h"
#include "particles.h"


void directdriver(struct particles *sources, struct particles *targets,
                  char *kernelName, double kernel_parameter, char *singularityHandling,
                  char *approximationName, double *pointwisePotential)
{
    int rank, numProcs, ierr;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    int numSources = sources->num;
    int numTargets = targets->num;

    double *source_x = sources->x;
    double *source_y = sources->y;
    double *source_z = sources->z;
    double *source_q = sources->q;
    double *source_w = sources->w;

    double *target_x = targets->x;
    double *target_y = targets->y;
    double *target_z = targets->z;
    double *target_q = targets->q;

    int numSourcesOnProc[numProcs];
    MPI_Allgather(&numSources, 1, MPI_INT, numSourcesOnProc, 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Win win_sources_x, win_sources_y, win_sources_z, win_sources_q, win_sources_w;
    MPI_Win_create(source_x, numSources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_x);
    MPI_Win_create(source_y, numSources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_y);
    MPI_Win_create(source_z, numSources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_z);
    MPI_Win_create(source_q, numSources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_q);
    MPI_Win_create(source_w, numSources*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win_sources_w);


    for (int procID = 1; procID < numProcs; ++procID) {

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

        Interaction_Direct_Compute(remote_sources->x, remote_sources->y, remote_sources->z,
                                   remote_sources->q, remote_sources->w,
                                   target_x, target_y, target_z, target_q,
                                   pointwisePotential, numSources, numTargets,
                                   kernelName, kernel_parameter, singularityHandling,
                                   approximationName);

        Particles_FreeSources(remote_sources);
    }

    //compute direct
    Interaction_Direct_Compute(source_x, source_y, source_z, source_q, source_w,
                               target_x, target_y, target_z, target_q,
                               pointwisePotential, numSources, numTargets,
                               kernelName, kernel_parameter, singularityHandling,
                               approximationName);

    //add correction
    Interaction_SubtractionPotentialCorrection(pointwisePotential, target_q, numTargets,
                               kernelName, kernel_parameter, singularityHandling);


    return;

} /* END of function pc_treecode */
