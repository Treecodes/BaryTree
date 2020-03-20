#ifndef H_BATCHFUNCTIONS_H
#define H_BATCHFUNCTIONS_H

#include "../particles/struct_particles.h"
#include "struct_tree.h"

void Batches_Alloc(struct Tree **batches, double *batch_lim,
                   struct Particles *particles, int batch_size);

void Batches_AllocArray(struct Tree **batches, int length);

void Batches_ReallocArray(struct Tree *batches, int length);

void Batches_Free(struct Tree *batches);

void Batches_Free_Win(struct Tree *batches);

void Batches_CreateTargetBatches(struct Tree *batches, struct Particles *particles,
                   int ibeg, int iend, int maxparnode, double *xyzmm);

void Batches_CreateSourceBatches(struct Tree *batches, struct Particles *particles,
                   int ibeg, int iend, int maxparnode, double *xyzmm);

#endif
