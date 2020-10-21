#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "../utilities/array.h"

#include "struct_particles.h"
#include "particles.h"


void Particles_Alloc(struct Particles **sources_addr, int length)
{
    *sources_addr = malloc(sizeof(struct Particles));
    struct Particles *sources = *sources_addr;

	sources->num = length;
    sources->x = NULL;
    sources->y = NULL;
    sources->z = NULL;
    sources->q = NULL;
    sources->w = NULL;
    
    if (sources->num > 0) {
        make_vector(sources->x, sources->num);
        make_vector(sources->y, sources->num);
        make_vector(sources->z, sources->num);
        make_vector(sources->q, sources->num);
        make_vector(sources->w, sources->num);
    }

    return;
}



void Particles_Free(struct Particles **sources_addr)
{
    struct Particles *sources = *sources_addr;

    if (sources != NULL) {
	    if (sources->x != NULL) free_vector(sources->x);
	    if (sources->y != NULL) free_vector(sources->y);
	    if (sources->z != NULL) free_vector(sources->z);
	    if (sources->q != NULL) free_vector(sources->q);
	    if (sources->w != NULL) free_vector(sources->w);
        free(sources);
    }
    
    sources = NULL;

    return;
}



void Particles_Targets_Reorder(struct Particles *targets, double *potential)
{
    int numpars = targets->num;
    int *reorder = targets->order;

    double *temp_energy;
    make_vector(temp_energy, numpars);
    for (int i = 0; i < numpars; i++) temp_energy[i] = potential[i];
    for (int i = 0; i < numpars; i++) potential[reorder[i]-1] = temp_energy[i];
    free_vector(temp_energy);

    double *temp_x;
    make_vector(temp_x, numpars);
    for (int i = 0; i < numpars; i++) temp_x[i] = targets->x[i];
    for (int i = 0; i < numpars; i++) targets->x[reorder[i]-1] = temp_x[i];
    free_vector(temp_x);

    double *temp_y;
    make_vector(temp_y, numpars);
    for (int i = 0; i < numpars; i++) temp_y[i] = targets->y[i];
    for (int i = 0; i < numpars; i++) targets->y[reorder[i]-1] = temp_y[i];
    free_vector(temp_y);

    double *temp_z;
    make_vector(temp_z, numpars);
    for (int i = 0; i < numpars; i++) temp_z[i] = targets->z[i];
    for (int i = 0; i < numpars; i++) targets->z[reorder[i]-1] = temp_z[i];
    free_vector(temp_z);

    double *temp_q;
    make_vector(temp_q, numpars);
    for (int i = 0; i < numpars; i++) temp_q[i] = targets->q[i];
    for (int i = 0; i < numpars; i++) targets->q[reorder[i]-1] = temp_q[i];
    free_vector(temp_q);

    return;
}



void Particles_Sources_Reorder(struct Particles *sources)
{
    int numpars = sources->num;
    int *reorder = sources->order;

    double *temp_x;
    make_vector(temp_x, numpars);
    for (int i = 0; i < numpars; i++) temp_x[i] = sources->x[i];
    for (int i = 0; i < numpars; i++) sources->x[reorder[i]-1] = temp_x[i];
    free_vector(temp_x);

    double *temp_y;
    make_vector(temp_y, numpars);
    for (int i = 0; i < numpars; i++) temp_y[i] = sources->y[i];
    for (int i = 0; i < numpars; i++) sources->y[reorder[i]-1] = temp_y[i];
    free_vector(temp_y);

    double *temp_z;
    make_vector(temp_z, numpars);
    for (int i = 0; i < numpars; i++) temp_z[i] = sources->z[i];
    for (int i = 0; i < numpars; i++) sources->z[reorder[i]-1] = temp_z[i];
    free_vector(temp_z);

    double *temp_q;
    make_vector(temp_q, numpars);
    for (int i = 0; i < numpars; i++) temp_q[i] = sources->q[i];
    for (int i = 0; i < numpars; i++) sources->q[reorder[i]-1] = temp_q[i];
    free_vector(temp_q);

    double *temp_w;
    make_vector(temp_w, numpars);
    for (int i = 0; i < numpars; i++) temp_w[i] = sources->w[i];
    for (int i = 0; i < numpars; i++) sources->w[reorder[i]-1] = temp_w[i];
    free_vector(temp_w);
    
    return;
}



void Particles_ConstructOrder(struct Particles *particles)
{
    make_vector(particles->order, particles->num);
    for (int i = 0; i < particles->num; i++) particles->order[i] = i+1;
    
    return;
}



void Particles_FreeOrder(struct Particles *particles)
{
    if (particles->order != NULL) free_vector(particles->order);
    
    return;
}



void Particles_Validate(struct Particles *sources, struct Particles *targets, struct RunParams *run_params)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    sources->sources_w_dummy = 0;
    targets->targets_q_dummy = 0;
        
    targets->targets_x_duplicate = 0;
    targets->targets_y_duplicate = 0;
    targets->targets_z_duplicate = 0;
    targets->targets_q_duplicate = 0;
    
    if (sources->x == NULL || sources->y == NULL || sources->z == NULL || sources->q == NULL
     || targets->x == NULL || targets->y == NULL || targets->z == NULL) {
        printf("[BaryTree]\n");
        printf("[BaryTree] ERROR! One or more required particle arrays is unset. Exiting.\n");
        printf("[BaryTree]\n");
        exit(1);
    }
    
    if (sources->w == NULL && run_params->singularity == SUBTRACTION) {
        printf("[BaryTree]\n");
        printf("[BaryTree] ERROR! Singularity subtraction requires weights for the sources. Exiting.\n");
        printf("[BaryTree]\n");
        exit(1);
    }
    
    if (sources->w == NULL) {
        make_vector(sources->w, sources->num);
        sources->sources_w_dummy = 1;
    }
    
    if (targets->q == NULL) {
        make_vector(targets->q, targets->num);
        targets->targets_q_dummy = 1;
    }
    
    if (sources->x == targets->x) {
        make_vector(targets->x, targets->num);
        memcpy(targets->x, sources->x, targets->num * sizeof(double));
        targets->targets_x_duplicate = 1;
    }

    if (sources->y == targets->y) {
        make_vector(targets->y, targets->num);
        memcpy(targets->y, sources->y, targets->num * sizeof(double));
        targets->targets_y_duplicate = 1;
    }

    if (sources->z == targets->z) {
        make_vector(targets->z, targets->num);
        memcpy(targets->z, sources->z, targets->num * sizeof(double));
        targets->targets_z_duplicate = 1;
    }

    if (sources->q == targets->q) {
        make_vector(targets->q, targets->num);
        memcpy(targets->q, sources->q, targets->num * sizeof(double));
        targets->targets_q_duplicate = 1;
    }
    
    if (targets->targets_x_duplicate || targets->targets_x_duplicate || targets->targets_z_duplicate
     || targets->targets_q_duplicate) {
        printf("[BaryTree]\n");
        printf("[BaryTree] Sources and targets cannot be the same location in memory.\n");
        printf("[BaryTree] Making duplicate arrays for targets as necessary.\n");
        printf("[BaryTree]\n");
    }

    return;
}
