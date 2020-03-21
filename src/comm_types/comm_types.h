#ifndef H_COMM_TYPES_FUNCTIONS_H
#define H_COMM_TYPES_FUNCTIONS_H

#include "../tree/struct_tree.h"
#include "../run_params/struct_run_params.h"

#include "struct_comm_types.h"


void CommTypesAndTrees_Construct(struct CommTypes **comm_types_addr, struct Tree ***let_trees_addr,
                                 struct Tree *tree, struct Tree *batches,
                                 struct RunParams *run_params);

void CommTypesAndTrees_Free(struct CommTypes **comm_types_addr, struct Tree ***let_trees_addr);


#endif /* H_COMM_TYPES_FUNCTIONS_H */
