#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>

#include "array.h"
#include "struct_kernel.h"

#include "kernel.h"

/* Routines to allocate, set, and free the kernel struct. */

void Kernel_Allocate(struct kernel *kernel, int numberOfParameters, char *name)
{

	kernel->name = name;
    kernel->numberOfParameters = numberOfParameters;
	if (kernel->numberOfParameters > 0) make_vector(kernel->parameters, numberOfParameters);

    return;
}



void Kernel_SetParams(struct kernel *kernel, double *parameters)
{
    for (int i = 0; i < kernel->numberOfParameters; i++) {
        kernel->parameters[i] = parameters[i];
    }

}



void Kernel_Free(struct kernel *kernel)
{
    if (kernel != NULL) {
	    if (kernel->numberOfParameters > 0) free_vector(kernel->parameters);
        free(kernel);
    }

    return;
}
