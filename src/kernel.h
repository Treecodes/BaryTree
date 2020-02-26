#ifndef H_KERNELFUNCTIONS_H
#define H_KERNELFUNCTIONS_H

#include "struct_kernel.h"


void Kernel_Allocate(struct kernel *kernel, int numberOfParameters, char *name);
void Kernel_SetParams(struct kernel *kernel, double *parameters);
void Kernel_Free(struct kernel *kernel);


#endif
