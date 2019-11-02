#include <math.h>
#include "kernels.h"
#include <float.h>


double coulombKernel( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa){

	double dx, dy, dz, r;

	dx = targetX-sourceX;
	dy = targetY-sourceY;
	dz = targetZ-sourceZ;

	r = sqrt( (dx*dx) + (dy*dy) + (dz*dz));


	if (r>DBL_MIN){
		return sourceQ*sourceW/r;
	}else{
		return 0.0;
	}
}

double yukawaKernel( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa){

	double dx, dy, dz, r;

	dx = targetX-sourceX;
	dy = targetY-sourceY;
	dz = targetZ-sourceZ;

	r = sqrt( (dx*dx) + (dy*dy) + (dz*dz));
	if (r>DBL_MIN){
		return sourceQ*sourceW*exp(-kappa*r)/r;
	}else{
		return 0.0;
	}
}

double coulombKernel_SS( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa){

	double dx, dy, dz, r, kappaSq;

	dx = targetX-sourceX;
	dy = targetY-sourceY;
	dz = targetZ-sourceZ;

	r = sqrt( (dx*dx) + (dy*dy) + (dz*dz));
	kappaSq = kappa*kappa;

	if (r>DBL_MIN){
		return (sourceQ - targetQ*exp(-r*r/kappaSq))*sourceW/r;
	}else{
		return 0.0;
	}
}

double yukawaKernel_SS( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa){

	double dx, dy, dz, r;

	dx = targetX-sourceX;
	dy = targetY-sourceY;
	dz = targetZ-sourceZ;

	r = sqrt( (dx*dx) + (dy*dy) + (dz*dz));


	if (r>DBL_MIN){
		return (sourceQ-targetQ)*sourceW*exp(-kappa*r)/r;
	}else{
		return 0.0;
	}
}
