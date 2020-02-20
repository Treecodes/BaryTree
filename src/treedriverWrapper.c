#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <limits.h>
#include <string.h>

#include "array.h"
#include "struct_kernel.h"
#include "kernel.h"
#include "struct_particles.h"
#include "particles.h"
#include "tools.h"
#include "tree.h"

#include "treedriver.h"
#include "treedriverWrapper.h"


/* definition of primary treecode driver */


void treedriverWrapper(int numTargets, int numSources,
		double *targetX, double *targetY, double *targetZ, double *targetValue,
		double *sourceX, double *sourceY, double *sourceZ, double *sourceValue, double *sourceWeight,
		double *outputArray, char *kernelName, int numberOfParameters, double * kernelParameters, char *singularityHandling, char *approximationName,
		int order, double theta, int maxparnode, int batch_size, int verbosity) {

    double sizeCheckFactor=1.0;
    if (strcmp(approximationName,"lagrange")==0){
        sizeCheckFactor=1.0;
    }else if (strcmp(approximationName,"hermite")==0){
        sizeCheckFactor=4.0;
    }

    struct kernel *kernel = NULL;
    kernel = malloc(sizeof (struct kernel));
    AllocateKernelStruct(kernel, numberOfParameters, kernelName);
    SetKernelParameters(kernel, kernelParameters);


	// Assemble the arrays of data into the particle structs.
	struct particles sources;
	struct particles targets;

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


	double time_tree[12];
	int tree_type = 1;   // particle cluster
	double tpeng = 0;

	// Initialize the potential
	for (int i=0; i<numTargets; i++){
            outputArray[i]=0.0;
    }



	// Call the treedriver
	treedriver(&sources, &targets,
			   order, theta, maxparnode, batch_size,
			   kernel, singularityHandling, approximationName, tree_type,
			   outputArray, time_tree, sizeCheckFactor, verbosity);
	MPI_Barrier(MPI_COMM_WORLD);


    FreeKernelStruct(kernel);

	return;
}
