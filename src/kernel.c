#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "array.h"
#include "struct_kernel.h"


/* Routines to allocate, set, and free the kernel struct. */

void AllocateKernelStruct(struct kernel *kernel, int numberOfParameters, char *name)
{

	kernel->name = name;
    kernel->numberOfParameters = numberOfParameters;
	make_vector(kernel->parameters, numberOfParameters);

    return;
}


void SetKernelParameters(struct kernel *kernel, double * parameters){

    for (int i=0;i<kernel->numberOfParameters;i++){
        kernel->parameters[i]=parameters[i];
    }

}



void FreeKernelStruct(struct kernel *kernel)
{
	free_vector(kernel->parameters);
    free(kernel);

    return;
}
