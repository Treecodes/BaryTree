#ifndef H_TREEDRIVERWRAPPER_H
#define H_TREEDRIVERWRAPPER_H

/* declaration of primary treecode driver wrapper */

void treedriverWrapper(int numTargets, int numSources,
		double *targetX, double *targetY, double *targetZ, double *targetValue,
		double *sourceX, double *sourceY, double *sourceZ, double *sourceValue, double *sourceWeight,
		double *outputArray, int pot_type, double kappa,
		int order, double theta, int maxparnode, int batch_size, MPI_Comm comm);

#endif /* H_TREEDRIVERWRAPPER_H */
