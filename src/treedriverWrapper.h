#ifndef H_TREEDRIVERWRAPPER_H
#define H_TREEDRIVERWRAPPER_H

/* declaration of primary treecode driver wrapper */

void treedriverWrapper(int numTargets, int numSources,
		double *targetX, double *targetY, double *targetZ, double *targetValue,
		double *sourceX, double *sourceY, double *sourceZ, double *sourceValue, double *sourceWeight,
		double *outputArray, char *kernelName, int numberOfParameters, double * kernelParameters, char *singularityHandling, char *approximationName,
		int order, double theta, int maxparnode, int batch_size, int verbosity);

#endif /* H_TREEDRIVERWRAPPER_H */
