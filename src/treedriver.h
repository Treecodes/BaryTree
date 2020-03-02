#ifndef H_TREEDRIVER_H
#define H_TREEDRIVER_H

#include "struct_particles.h"
#include "struct_run_params.h"


void treedriver(struct particles *sources, struct particles *targets, struct RunParams *run_params,
                double *potential_array, double *time_tree);


#endif /* H_TREEDRIVER_H */
