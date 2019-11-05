/* Interaction Kernels */
#ifndef H_KERNELS_H
#define H_KERNELS_H

double coulombKernel( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa);

double yukawaKernel( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa);

double coulombKernel_SS_direct( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa);

double coulombKernel_SS_approx( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa);

double yukawaKernel_SS_direct( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa);

double yukawaKernel_SS_approx( double targetX, double targetY, double targetZ, double targetQ,
					double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
					double kappa);


#endif /* H_KERNELS_H */
