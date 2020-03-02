#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#include "array.h"
#include "struct_kernel.h"
#include "kernel.h"
#include "struct_particles.h"
#include "struct_output.h"
#include "particles.h"
#include "tools.h"
#include "tree.h"

#include "treedriver.h"
#include "BaryTreeInterface.h"


/* definition of primary treecode driver */


void BaryTreeInterface(int numTargets, int numSources,
		double *targetX, double *targetY, double *targetZ, double *targetValue,
		double *sourceX, double *sourceY, double *sourceZ, double *sourceValue, double *sourceWeight,
		double *outputArray, char *kernelName, int numberOfParameters, double *kernelParameters,
        char *singularityHandling, char *approximationName,
		int interpOrder, double theta, int maxPerLeaf, int maxPerBatch, int verbosity)
{

	double timing[12];
	int treeType = 1;   // particle-cluster
    memset(outputArray, 0, numTargets * sizeof(double));

    double sizeCheckFactor = 1.0;
    if (strcmp(approximationName, "lagrange") == 0) {
        sizeCheckFactor = 1.0;
    } else if (strcmp(approximationName, "hermite") == 0) {
        sizeCheckFactor = 4.0;
    }

    struct kernel *kernel = malloc(sizeof (struct kernel));
    Kernel_Allocate(kernel, numberOfParameters, kernelName);
    Kernel_SetParams(kernel, kernelParameters);


	struct particles sources, targets;

	targets.num = numTargets;
	targets.x = targetX;
	targets.y = targetY;
	targets.z = targetZ;
	targets.q = targetValue;

	sources.num = numSources;
	sources.x = sourceX;
	sources.y = sourceY;
	sources.z = sourceZ;
	sources.q = sourceValue;
	sources.w = sourceWeight;

	struct output output;
	int forces=0;
	Output_Alloc(output, numTargets, forces);



	treedriver(&sources, &targets,
			   interpOrder, theta, maxPerLeaf, maxPerBatch,
			   kernel, singularityHandling, approximationName, treeType,
			   output, timing, sizeCheckFactor, verbosity);

	MPI_Barrier(MPI_COMM_WORLD);


    Kernel_Free(kernel);

	return;
}
