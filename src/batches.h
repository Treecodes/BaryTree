#ifndef H_BATCHFUNCTIONS_H
#define H_BATCHFUNCTIONS_H

#include "nodes_struct.h"
#include "particles_struct.h"

void Batches_Setup(struct tnode_array **batches, double *batch_lim,
                   struct particles *particles, int batch_size);

void Batches_CreateTargetBatches(struct tnode_array *batches, struct particles *particles,
                   int ibeg, int iend, int maxparnode, double *xyzmm);

void Batches_CreateSourceBatches(struct tnode_array *batches, struct particles *particles,
                   int ibeg, int iend, int maxparnode, double *xyzmm);

void reorder_targets_and_potential(struct particles *targets, double *tEn, int numpars);

#endif
