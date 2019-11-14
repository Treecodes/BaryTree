#include <math.h>
#include <float.h>

#include "kernels_direct.h"


#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double coulombKernel(double targetX, double targetY, double targetZ, double targetQ,
                     double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                     double kappa)
{
    double dx = targetX-sourceX;
    double dy = targetY-sourceY;
    double dz = targetZ-sourceZ;

    double r = sqrt(dx*dx + dy*dy + dz*dz);

    if (r > DBL_MIN)
        return sourceQ * sourceW / r;
    else
        return 0.0;
}




#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double yukawaKernel(double targetX, double targetY, double targetZ, double targetQ,
                    double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                    double kappa)
{
    double dx = targetX - sourceX;
    double dy = targetY - sourceY;
    double dz = targetZ - sourceZ;

    double r = sqrt(dx*dx + dy*dy + dz*dz);

    if (r > DBL_MIN)
        return sourceQ * sourceW * exp(-kappa*r) / r;
    else
        return 0.0;
}




#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double coulombKernel_SS(double targetX, double targetY, double targetZ, double targetQ,
                        double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                        double kappa)
{
    double dx = targetX - sourceX;
    double dy = targetY - sourceY;
    double dz = targetZ - sourceZ;

    double r = sqrt(dx*dx + dy*dy + dz*dz);
    double kappaSq = kappa*kappa;

    if (r > DBL_MIN)
        return (sourceQ - targetQ * exp(-r*r/kappaSq)) * sourceW / r;
    else
        return 0.0;
}




#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double yukawaKernel_SS(double targetX, double targetY, double targetZ, double targetQ,
                       double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                       double kappa)
{
    double dx = targetX - sourceX;
    double dy = targetY - sourceY;
    double dz = targetZ - sourceZ;

    double r = sqrt(dx*dx + dy*dy + dz*dz);
    double G = exp(-kappa*r) / r;

    if (r > DBL_MIN)
        return (sourceQ - targetQ) * sourceW * G;
    else
        return 0.0;
}
