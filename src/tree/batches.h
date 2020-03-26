#ifndef H_BATCH_FUNCTIONS_H
#define H_BATCH_FUNCTIONS_H

#include "../particles/struct_particles.h"
#include "../run_params/struct_run_params.h"

#include "struct_tree.h"


void Batches_Sources_Construct(struct Tree **batches_addr, struct Particles *sources,
                struct RunParams *run_params);
                
void Batches_Targets_Construct(struct Tree **batches_addr, struct Particles *targets,
                struct RunParams *run_params);

void Batches_Alloc(struct Tree **batches_addr, int length);

void Batches_Free(struct Tree **batches_addr);

void Batches_Free_Win(struct Tree **batches_addr);

void Batches_Print(struct Tree *batches);


#endif
