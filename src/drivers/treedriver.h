#ifndef H_TREEDRIVER_H
#define H_TREEDRIVER_H

#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"


void treedriver(struct Particles *sources, struct Particles *targets, struct RunParams *run_params,
                double *potential_array, double *time_tree);


#endif /* H_TREEDRIVER_H */
