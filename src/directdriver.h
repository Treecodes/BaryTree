#ifndef H_DIRECTDRIVER_H
#define H_DIRECTDRIVER_H

#include "struct_particles.h"
#include "struct_run_params.h"


void directdriver(struct particles *sources, struct particles *targets, struct RunParams *run_params,
                  double *potential_array, double *time_direct);


#endif /* H_DIRECTDRIVER_H */
