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


//	int particleOrder[numTargets];
//	for (int i=0; i<numTargets; i++){ particleOrder[i]=i;}  // should order start at 0 or 1?  Looks like 0, as in main.c

	// Assemble the arrays of data into the particle structs.
	struct particles *sources = NULL;
	struct particles *targets = NULL;
	sources = malloc(sizeof(struct particles));
	targets = malloc(sizeof(struct particles));

	targets->num = numTargets;
	targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));
    targets->w = malloc(targets->num*sizeof(double));  // need to allocated targets->w so that the free_particles function doesn't segfault.
    memcpy(targets->x, targetX, targets->num*sizeof(double));
    memcpy(targets->y, targetY, targets->num*sizeof(double));
    memcpy(targets->z, targetZ, targets->num*sizeof(double));
    memcpy(targets->q, targetValue, targets->num*sizeof(double));
//	targets->x = targetX;
//	targets->y = targetY;
//	targets->z = targetZ;
//	targets->q = targetValue;
//	targets->order = particleOrder;

	sources->num = numSources;
	sources->x = malloc(sources->num*sizeof(double));
    sources->y = malloc(sources->num*sizeof(double));
    sources->z = malloc(sources->num*sizeof(double));
    sources->q = malloc(sources->num*sizeof(double));
    sources->w = malloc(sources->num*sizeof(double));
    memcpy(sources->x, sourceX, sources->num*sizeof(double));
    memcpy(sources->y, sourceY, sources->num*sizeof(double));
    memcpy(sources->z, sourceZ, sources->num*sizeof(double));
    memcpy(sources->q, sourceValue, sources->num*sizeof(double));
    memcpy(sources->w, sourceWeight, sources->num*sizeof(double));
//	sources->x = sourceX;
//	sources->y = sourceY;
//	sources->z = sourceZ;
//	sources->q = sourceValue;
//	sources->w = sourceWeight;
//	sources->order = particleOrder;

	double time_tree[9];
	int tree_type = 1;   // particle cluster
	double tpeng = 0;

	// Initialize the potential

	for (int i=0; i<numTargets; i++){
            outputArray[i]=0.0;
    }



	// Call the treedriver
//    printf("In wrapper, sources->x exist before call to treedriver? %1.3e\n", sources->x[3]);
//    printf("In wrapper, sourceX exist before call to treedriver? %1.3e\n", sourceX[3]);
    MPI_Barrier(MPI_COMM_WORLD);
	treedriver(sources, targets,
			   order, theta, maxparnode, batch_size,
			   kernel, singularityHandling, approximationName, tree_type,
			   outputArray, time_tree, sizeCheckFactor, verbosity);
	MPI_Barrier(MPI_COMM_WORLD);


//    tpeng = sum(outputArray, targets->num); // this isn't being used, no need to sum the computed potential

    // free the particle structs (but not the member arrays themselves, which already existed before the call to treedriverwrapper and need to persist)

//	printf("In wrapper, sources->x still exist after call to treedriver? %1.3e\n", sources->x[3]);
//    printf("In wrapper, sourceX still exist after call to treedriver? %1.3e\n", sourceX[3]);

//    MPI_Barrier(MPI_COMM_WORLD);
//    free(sources);
//	free(targets);
//	MPI_Barrier(MPI_COMM_WORLD);
//    printf("In wrapper, sourceX still exist after call to free? %1.3e\n", sourceX[3]);  // sourceX still exists.  sources->x no longer exists.
//    MPI_Barrier(MPI_COMM_WORLD);


    Particles_FreeSources(sources);
    Particles_FreeSources(targets);
    FreeKernelStruct(kernel);

	return;
}
