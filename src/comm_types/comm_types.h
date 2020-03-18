#ifndef H_COMM_TYPES_FUNCTIONS_H
#define H_COMM_TYPES_FUNCTIONS_H

#include "../tree/struct_nodes.h"
#include "../run_params/struct_run_params.h"

#include "struct_comm_types.h"


void CommTypesAndTrees_Construct(struct CommTypes **comm_types_addr, struct tnode_array ***let_tree_arrays_addr,
                                 struct tnode_array *tree_array, struct tnode_array *batches,
                                 struct RunParams *run_params);

void CommTypesAndTrees_Free(struct CommTypes *comm_types, struct tnode_array **let_tree_arrays);


#endif /* H_COMM_TYPES_FUNCTIONS_H */
