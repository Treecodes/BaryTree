#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <limits.h>
#include <string.h>

#include "array.h"
#include "struct_particles.h"
#include "tools.h"
#include "tree.h"

#include "treedriver.h"
#include "treedriverWrapper.h"


/* definition of primary treecode driver */


void treedriverWrapper(int numTargets, int numSources,
		double *targetX, double *targetY, double *targetZ, double *targetValue,
		double *sourceX, double *sourceY, double *sourceZ, double *sourceValue, double *sourceWeight,
		double *outputArray, char *kernelName, double kappa, char *singularityHandling, char *approximationName,
		int order, double theta, int maxparnode, int batch_size) {

    int verbosity=0;
    double sizeCheckFactor=1.0;
//	int particleOrder[numTargets];
//	for (int i=0; i<numTargets; i++){ particleOrder[i]=i;}  // should order start at 0 or 1?  Looks like 0, as in main.c

	// Assemble the arrays of data into the particle structs.
	struct particles *sources = NULL;
	struct particles *targets = NULL;
	sources = malloc(sizeof(struct particles));
	targets = malloc(sizeof(struct particles));

	targets->num = numTargets;
	targets->x = targetX;
	targets->y = targetY;
	targets->z = targetZ;
	targets->q = targetValue;
//	targets->order = particleOrder;

	sources->num = numSources;
	sources->x = sourceX;
	sources->y = sourceY;
	sources->z = sourceZ;
	sources->q = sourceValue;
	sources->w = sourceWeight;
//	sources->order = particleOrder;

	double time_tree[9];
	int tree_type = 1;   // particle cluster
	double tpeng = 0;

	// Initialize the potential
	if (strcmp(singularityHandling,"skipping")==0){
	    for (int i=0; i<numTargets; i++){
            outputArray[i]=0.0;
        }
	}else if (strcmp(singularityHandling,"subtraction")==0){

	    if (strcmp(kernelName,"coulomb")==0){
	        if (verbosity>0) printf("Initializing for coulomb singularity subtraction.\n");
	        for (int i=0; i<numTargets; i++){
                outputArray[i]=0.0;
	        }
        } else if (strcmp(kernelName,"yukawa")==0){
            if (verbosity>0) printf("Initializing for yukawa singularity subtraction.\n");
            for (int i=0; i<numTargets; i++){
                outputArray[i]=0.0;
            }
	    } else{
	        printf("Not sure how to initialize outputArray.  What is the kernel?\n");
	        exit(-1);
	    }
	} else {
	    printf("Not sure how to initialize outputArray.  How is singularity being handled?\n");
        exit(-1);
	}


	// Call the treedriver
	treedriver(sources, targets,
			   order, theta, maxparnode, batch_size,
			   kernelName, kappa, singularityHandling, approximationName, tree_type,
			   outputArray, time_tree, sizeCheckFactor, verbosity);

//    tpeng = sum(outputArray, targets->num); // this isn't being used, no need to sum the computed potential

    // free the particle structs (but not the member arrays themselves, which already existed before the call to treedriverwrapper and need to persist)
	free(sources);
	free(targets);

	return;

}
