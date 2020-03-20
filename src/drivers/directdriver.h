#ifndef H_DIRECTDRIVER_H
#define H_DIRECTDRIVER_H

#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"


void directdriver(struct Particles *sources, struct Particles *targets, struct RunParams *run_params,
                  double *potential_array, double *time_direct);


#endif /* H_DIRECTDRIVER_H */
