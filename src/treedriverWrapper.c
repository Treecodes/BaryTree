#include <stdio.h>
#include <mpi.h>
#include <limits.h>

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
		double *outputArray, int pot_type, double kappa,
		int order, double theta, int maxparnode, int batch_size, int numDevices, int numThreads) {

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
	targets->order=particleOrder;


	sources->num = numSources;
	sources->x = sourceX;
	sources->y = sourceY;
	sources->z = sourceZ;
	sources->q = sourceValue;
	sources->w = sourceWeight;
	sources->order=particleOrder;

//	for (int j=0;j<15;j++){printf("target->order[%d] = %d\n", j, targets->order[j]);}


	// Set other things, like *timetree and pot_type...

	double time_tree[4];
	int tree_type=1;   // particle cluster
	double tpeng = 0;

	// Initialize all GPUs
	if (numDevices>0){
		#pragma omp parallel num_threads(numDevices)
			{
			acc_set_device_num(omp_get_thread_num(),acc_get_device_type());
			acc_init(acc_get_device_type());
			}
	}

	// Call the treedriver
	treedriver(sources, targets,
			   order, theta, maxparnode, batch_size,
			   pot_type, kappa, tree_type,
			   outputArray, &tpeng, time_tree, numDevices, numThreads);


//	free_vector(sources->x);
//	free_vector(sources->y);
//	free_vector(sources->z);
//	free_vector(sources->q);
//	free_vector(sources->w);
//	free_vector(sources->order);
	free(sources);

//	free_vector(targets->x);
//	free_vector(targets->y);
//	free_vector(targets->z);
//	free_vector(targets->q);
//	free_vector(targets->w);
//	free_vector(targets->order);
	free(targets);

	return;

}
