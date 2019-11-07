#include <stdio.h>
#include <mpi.h>
#include <limits.h>
#include<string.h>

#include "array.h"
#include "particles.h"
#include "tools.h"
#include "tree.h"

#include "treedriver.h"
#include "treedriverWrapper.h"


/* definition of primary treecode driver */


void treedriverWrapper(int numTargets, int numSources,
		double *targetX, double *targetY, double *targetZ, double *targetValue,
		double *sourceX, double *sourceY, double *sourceZ, double *sourceValue, double *sourceWeight,
		double *outputArray, char *kernelName, double kappa,
		int order, double theta, int maxparnode, int batch_size) {

	// Set up kernels
	if       (strcmp(kernelName,"coulomb")==0){

	}else if (strcmp(kernelName,"yukawa")==0){

	}else if (strcmp(kernelName,"coulomb_SS")==0){

	}else if (strcmp(kernelName,"yukawa_SS")==0){

	}else{
		return;
	}


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

	// Call the treedriver
	treedriver(sources, targets,
			   order, theta, maxparnode, batch_size,
			   kernelName, kappa, tree_type,
			   outputArray, &tpeng, time_tree);

	free(sources);
	free(targets);

	return;

}
