#ifndef H_KERNELFUNCTIONS_H
#define H_KERNELFUNCTIONS_H

#include "struct_kernel.h"


void AllocateKernelStruct(struct kernel *kernel, int numberOfParameters, char *name);
void SetKernelParameters(struct kernel *kernel, double * parameters);
void FreeKernelStruct(struct kernel *kernel);


#endif
