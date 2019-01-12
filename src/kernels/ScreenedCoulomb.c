#include <math.h>

void screenedCoulombSingularitySubtractionKernel( int numberOfTargets, int numberOfInterpolationPoints, int indexOfFirstTarget,
					double *targetX, double *targetY, double *targetZ, double *targetVal,
					double *interpolationX, double interpolationY, double *interpolationZ, double *interpolationVal,
					double screeningK, double *kernelMatrix){

	// indexOfFirstTarget = batch_ind[0] - 1 with current convention

	int i,j;
	double tx, ty, tz;
	double dx, dy, dz, r;

	for (i = 0; i < numberOfTargets; i++){
		tx = targetX[ indexOfFirstTarget + i];
		ty = targetY[ indexOfFirstTarget + i];
		tz = targetZ[ indexOfFirstTarget + i];

		for (j = 0; j < numberOfInterpolationPoints; j++){

			// Compute x, y, and z distances between target i and interpolation point j
			dx = targetX[ indexOfFirstTarget + i] - interpolationX[j];
			dy = targetY[ indexOfFirstTarget + i] - interpolationY[j];
			dz = targetZ[ indexOfFirstTarget + i] - interpolationZ[j];

			// Evaluate Kernel, store in kernelMatrix[i][j]
			r = sqrt( dx*dx + dy*dy + dz*dz);
			kernelMatrix[i*numberOfInterpolationPoints + j] = interpolationVal[j] * exp(-screeningK*r) / r;

		}

	}

	return;
}
