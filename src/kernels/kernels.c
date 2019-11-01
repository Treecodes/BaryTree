#include <math.h>
#include "kernels.h"

double coulombKernel( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa){

	double dx, dy, dz, r;

	dx = targetX-sourceX;
	dy = targetY-sourceY;
	dz = targetZ-sourceZ;

	r = sqrt( (dx*dx) + (dy*dy) + (dz*dz));

	return sourceQ*sourceW/r;
}

double yukawaKernel( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa){

	double dx, dy, dz, r;

	dx = targetX-sourceX;
	dy = targetY-sourceY;
	dz = targetZ-sourceZ;

	r = sqrt( (dx*dx) + (dy*dy) + (dz*dz));

	return sourceQ*sourceW*exp(-kappa*r)/r;
}
