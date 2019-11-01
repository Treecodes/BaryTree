#include <math.h>

double coulombKernel( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW){

	double dx, dy, dz, r;

	dx = targetX-sourceX;
	dy = targetY-sourceY;
	dz = targetZ-sourceZ;

	r = sqrt( (dx*dx) + (dy*dy) + (dz*dz));

//	*KernelValue = sourceQ*sourceW/r;

	return sourceQ*sourceW/r;
}
