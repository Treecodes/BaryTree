#ifndef H_PARTICLEFUNCTIONS_H
#define H_PARTICLEFUNCTIONS_H

#include "struct_particles.h"

void Particles_Alloc(struct particles **sources_addr, int length);

void Particles_Free(struct particles *sources);

void Particles_ReorderTargetsAndPotential(struct particles *targets, double *tEn);

#endif
