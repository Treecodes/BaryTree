#ifndef H_PARTICLEFUNCTIONS_H
#define H_PARTICLEFUNCTIONS_H

#include "struct_particles.h"

void Particles_Alloc(struct Particles **sources_addr, int length);

void Particles_Free(struct Particles *sources);

void Particles_ReorderTargetsAndPotential(struct Particles *targets, double *tEn);

#endif
