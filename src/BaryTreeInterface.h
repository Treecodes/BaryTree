#ifndef H_BARYTREE_INTERFACE_H
#define H_BARYTREE_INTERFACE_H

void BaryTreeInterface(int numTargets, int numSources,
		double *targetX, double *targetY, double *targetZ, double *targetValue,
		double *sourceX, double *sourceY, double *sourceZ, double *sourceValue, double *sourceWeight,
		double *outputArray, char *kernelName, int numberOfParameters, double *kernelParameters,
        char *singularityHandling, char *approximationName,
		int interpOrder, double theta, int maxPerLeaf, int maxPerBatch, int verbosity);

#endif /* H_BARYTREE_INTERFACE_H */
