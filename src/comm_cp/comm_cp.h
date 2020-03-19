#ifndef H_COMM_CP_FUNCTIONS_H
#define H_COMM_CP_FUNCTIONS_H

#include "../tree/struct_nodes.h"
#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"


void Comm_CP_ConstructAndGetData(struct tnode_array **remote_batches_addr, struct particles **remote_sources_addr,
                                 const struct tnode_array *tree_array, const struct tnode_array *batches,
                                 const struct particles *sources, const struct RunParams *run_params);


#endif /* H_COMM_CP_FUNCTIONS_H */
