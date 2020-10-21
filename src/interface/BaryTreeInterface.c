#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <limits.h>
#include <string.h>


#include "../particles/struct_particles.h"
#include "../particles/particles.h"

#include "../run_params/struct_run_params.h"
#include "../run_params/run_params.h"

#include "../drivers/treedriver.h"
#include "BaryTreeInterface.h"


void BaryTreeInterface(int numTargets, int numSources,
		double *targetX, double *targetY, double *targetZ, double *targetValue,
		double *sourceX, double *sourceY, double *sourceZ, double *sourceValue, double *sourceWeight,
		double *outputArray,
        KERNEL kernel, int numKernelParams, double *kernelParams,
        SINGULARITY singularity, APPROXIMATION approximation, COMPUTE_TYPE compute_type,
		double theta, int interpDegree, int maxPerSourceLeaf, int maxPerTargetLeaf,
        double sizeCheck, double beta, int verbosity)
{

	double timing[12];
    memset(outputArray, 0, numTargets * sizeof(double));

    struct RunParams *run_params = NULL;
    RunParams_Setup(&run_params,
                    kernel, numKernelParams, kernelParams,
                    approximation, singularity, compute_type,
                    theta, interpDegree,
                    maxPerSourceLeaf, maxPerTargetLeaf, sizeCheck,
                    beta, verbosity);

	struct Particles sources, targets;

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


	treedriver(&sources, &targets, run_params, outputArray, timing);
	MPI_Barrier(MPI_COMM_WORLD);


    RunParams_Free(&run_params);

	return;
}
