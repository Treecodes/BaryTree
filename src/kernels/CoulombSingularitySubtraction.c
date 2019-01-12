#include <math.h>


void coulombSingularitySubtractionKernel( int numberOfTargets, int numberOfInterpolationPoints, int indexOfFirstTarget,
					double *targetX, double *targetY, double *targetZ, double *targetVal,
					double *interpolationX, double interpolationY, double *interpolationZ, double *interpolationVal,
					double GaussianSmoothingAlphaSq, double *kernelMatrix){

	// indexOfFirstTarget = batch_ind[0] - 1 with current convention

	int i,j;
	double tx, ty, tz, targetValue;
	double dx, dy, dz, r;

	for (i = 0; i < numberOfTargets; i++){
		tx = targetX[ indexOfFirstTarget + i];
		ty = targetY[ indexOfFirstTarget + i];
		tz = targetZ[ indexOfFirstTarget + i];
		targetValue = targetVal[indexOfFirstTarget + i];

		for (j = 0; j < numberOfInterpolationPoints; j++){

			// Compute x, y, and z distances between target i and interpolation point j
			dx = targetX[ indexOfFirstTarget + i] - interpolationX[j];
			dy = targetY[ indexOfFirstTarget + i] - interpolationY[j];
			dz = targetZ[ indexOfFirstTarget + i] - interpolationZ[j];
			r = sqrt( dx*dx + dy*dy + dz*dz);

			// Evaluate Kernel, store in kernelMatrix[i][j]
			kernelMatrix[i*numberOfInterpolationPoints + j] = ( interpolationVal[j] - targetValue* exp(- r*r/GaussianSmoothingAlphaSq) ) / r;
            // Python version: V_Hartree_new[globalID] += weight_s * (rho_s -   rho_t * exp(- r*r / alphasq )   ) / r  # increment the new wavefunction value

		}

	}

	return;
}
