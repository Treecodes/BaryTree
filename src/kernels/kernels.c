#include <math.h>
#include "kernels.h"
#include <float.h>


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
double coulombKernel_SS_direct(double targetX, double targetY, double targetZ, double targetQ,
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
double coulombKernel_SS_approx(double targetX, double targetY, double targetZ, double targetQ,
                               double clusterX, double clusterY, double clusterZ, double clusterQ, double clusterW,
                               double kappa)
{
    double dx = targetX - clusterX;
    double dy = targetY - clusterY;
    double dz = targetZ - clusterZ;

    double r = sqrt(dx*dx + dy*dy + dz*dz);
    double kappaSq = kappa*kappa;

    if (r > DBL_MIN)
        return (clusterQ - targetQ * clusterW * exp(-r*r/kappaSq)) / r;
    else
        return 0.0;
}


#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double yukawaKernel_SS_direct(double targetX, double targetY, double targetZ, double targetQ,
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


#ifdef OPENACC_ENABLED
#pragma acc routine seq
#endif
double yukawaKernel_SS_approx(double targetX, double targetY, double targetZ, double targetQ,
                              double sourceX, double sourceY, double sourceZ, double sourceQ, double sourceW,
                              double kappa)
{
    double dx = targetX - sourceX;
    double dy = targetY - sourceY;
    double dz = targetZ - sourceZ;

    double r = sqrt(dx*dx + dy*dy + dz*dz);
    double G = exp(-kappa*r) / r;

    if (r > DBL_MIN)
        return (sourceQ - targetQ * sourceW) * G;
    else
        return 0.0;
}
