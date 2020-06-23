#ifndef H_BARYTREE_INTERFACE_H
#define H_BARYTREE_INTERFACE_H

#include "../utilities/enums.h"


void BaryTreeInterface(int numTargets, int numSources,
		double *targetX, double *targetY, double *targetZ, double *targetValue,
		double *sourceX, double *sourceY, double *sourceZ, double *sourceValue, double *sourceWeight,
		double *outputArray,
		double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
        KERNEL kernel, int numKernelParams, double *kernelParams,
        SINGULARITY singularity, APPROXIMATION approximation, COMPUTE_TYPE compute_type,
		int interpOrder, double theta, int maxPerSourceLeaf, int maxPerTargetLeaf,
        double sizeCheck, int verbosity);


#endif /* H_BARYTREE_INTERFACE_H */
