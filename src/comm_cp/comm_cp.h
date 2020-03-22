#ifndef H_COMM_CP_FUNCTIONS_H
#define H_COMM_CP_FUNCTIONS_H

#include "../tree/struct_tree.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"


void Comm_CP_ConstructAndGetData(struct Tree **remote_batches_addr, struct Particles **remote_sources_addr,
                                 const struct Tree *tree_array, const struct Tree *batches,
                                 const struct Particles *sources, const struct RunParams *run_params);


#endif /* H_COMM_CP_FUNCTIONS_H */
