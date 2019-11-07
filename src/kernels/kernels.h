/* Interaction Kernels */
#ifndef H_KERNELS_H
#define H_KERNELS_H
 
#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double coulombKernel(double targetX, double targetY, double targetZ, double targetQ,
                    double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                    double kappa);

#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double yukawaKernel(double targetX, double targetY, double targetZ, double targetQ,
                    double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                    double kappa);

#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double coulombKernel_SS_direct(double targetX, double targetY, double targetZ, double targetQ,
                    double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                    double kappa);

#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double coulombKernel_SS_approx(double targetX, double targetY, double targetZ, double targetQ,
                    double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                    double kappa);

#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double yukawaKernel_SS_direct(double targetX, double targetY, double targetZ, double targetQ,
                    double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                    double kappa);

#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double yukawaKernel_SS_approx(double targetX, double targetY, double targetZ, double targetQ,
                    double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                    double kappa);


#endif /* H_KERNELS_H */
