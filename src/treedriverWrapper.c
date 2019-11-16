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


	int particleOrder[numTargets];
	for (int i=0; i<numTargets; i++){ particleOrder[i]=i;}  // should order start at 0 or 1?  Looks like 0, as in main.c

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
	targets->order = particleOrder;

	sources->num = numSources;
	sources->x = sourceX;
	sources->y = sourceY;
	sources->z = sourceZ;
	sources->q = sourceValue;
	sources->w = sourceWeight;
	sources->order = particleOrder;

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
//	        printf("Initializing for coulomb singularity subtraction.\n");
	        for (int i=0; i<numTargets; i++){
                outputArray[i]=0.0;
//                outputArray[i]=2.0*M_PI*kappa*kappa*targets->q[i];
	        }
        } else if (strcmp(kernelName,"yukawa")==0){
//            printf("Initializing for yukawa singularity subtraction.\n");
            for (int i=0; i<numTargets; i++){
                outputArray[i]=0.0;
//                outputArray[i]=4.0*M_PI*targets->q[i]/kappa/kappa;  // 4*pi*f_t/k**2
            }
	    } else{
	        printf("Not sure how to initialize outputArray.  What is the kernel?\n");
	        exit(1);
	    }
	} else {
	    printf("Not sure how to initialize outputArray.  How is singularity being handled?\n");
        exit(1);
	}


	// Call the treedriver
	treedriver(sources, targets,
			   order, theta, maxparnode, batch_size,
			   kernelName, kappa, singularityHandling, approximationName, tree_type,
			   outputArray, time_tree);

    tpeng = sum(outputArray, targets->num);

	free(sources);
	free(targets);

	return;

}
