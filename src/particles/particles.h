#ifndef H_PARTICLE_FUNCTIONS_H
#define H_PARTICLE_FUNCTIONS_H

#include "../run_params/struct_run_params.h"
#include "struct_particles.h"

void Particles_Alloc(struct Particles **particles_addr, int length);

void Particles_Free(struct Particles **particles_addr);

void Particles_Targets_Reorder(struct Particles *targets, double *potential);

void Particles_Sources_Reorder(struct Particles *sources);

void Particles_ConstructOrder(struct Particles *particles);

void Particles_FreeOrder(struct Particles *particles);

void Particles_Validate(struct Particles *sources, struct Particles *targets,
        struct RunParams *run_params);


#endif /* H_PARTICLE_FUNCTIONS */
