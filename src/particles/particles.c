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



void Particles_Free(struct Particles *sources)
{
    if (sources != NULL) {
	    if (sources->x != NULL) free_vector(sources->x);
	    if (sources->y != NULL) free_vector(sources->y);
	    if (sources->z != NULL) free_vector(sources->z);
	    if (sources->q != NULL) free_vector(sources->q);
	    if (sources->w != NULL) free_vector(sources->w);
        free(sources);
    }

    return;
}



void Particles_ReorderTargetsAndPotential(struct Particles *targets, double *tEn)
{
    int numpars = targets->num;
    int *reorder = targets->order;

    double *temp_energy;
    make_vector(temp_energy, numpars);
    for (int i = 0; i < numpars; i++) temp_energy[i] = tEn[i];
    for (int i = 0; i < numpars; i++) tEn[reorder[i]-1] = temp_energy[i];
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

}
