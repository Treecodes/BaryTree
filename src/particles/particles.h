#ifndef H_PARTICLEFUNCTIONS_H
#define H_PARTICLEFUNCTIONS_H

#include "struct_particles.h"

void Particles_AllocSources(struct particles *sources, int length);

void Particles_FreeSources(struct particles *sources);

void Particles_ReorderTargetsAndPotential(struct particles *targets, double *tEn);

#endif
