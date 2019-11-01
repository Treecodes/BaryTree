#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <float.h>

#include "array.h"
#include "tools.h"
#include "direct.h"
#include "kernels/kernels.h"


void direct_eng(double *xS, double *yS, double *zS, double *qS, double *wS,
				double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa)
{
    /* local variables */
    int i, j;
    double tx, ty, tz, xi, yi, zi, qi, teng, rad;

#ifdef OPENACC_ENABLED
    #pragma acc data copyin (xS[0:numparsS], yS[0:numparsS], zS[0:numparsS], \
                             qS[0:numparsS], wS[0:numparsS], xT[0:numparsT], \
                             yT[0:numparsT], zT[0:numparsT], qT[0:numparsT])
    {
#endif

    if (pot_type == 0 || pot_type == 4) { // Coulomb with singularity skipping.  Lagrange or Hermite.
#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (i = 0; i < numparsT; i++) {
            xi = xT[i];
            yi = yT[i];
            zi = zT[i];
            teng = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (j = 0; j < numparsS; j++) {
                tx = xi - xS[j];
                ty = yi - yS[j];
                tz = zi - zS[j];
                rad = sqrt(tx*tx + ty*ty + tz*tz);
                if (rad > 1e-14) {
//                    teng = teng + qS[j] * wS[j] / rad;
                	teng = teng + coulombKernel(xi, yi, zi, 0.0, xS[j], yS[j], zS[j], qS[j], wS[j]);
                }
            }
            denergy[i] +=  teng;
        }
#ifdef OPENACC_ENABLED
        }
#endif

    } else if (pot_type == 1 || pot_type == 5) {
#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (i = 0; i < numparsT; i++) {
            xi = xT[i];
            yi = yT[i];
            zi = zT[i];
            teng = 0.0;

#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (j = 0; j < numparsS; j++) {
                tx = xi - xS[j];
                ty = yi - yS[j];
                tz = zi - zS[j];
                rad = sqrt(tx*tx + ty*ty + tz*tz);
                if (rad > 1e-14) {
                    teng = teng + qS[j] * wS[j] * exp(-kappa * rad) / rad;
                }
            }
            denergy[i] += teng;
        }
#ifdef OPENACC_ENABLED
        }
#endif

    } else if (pot_type == 2 || pot_type == 6) {
        double kappaSq = kappa*kappa;
#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (i = 0; i < numparsT; i++) {
            xi = xT[i];
            yi = yT[i];
            zi = zT[i];
            qi = qT[i];
            teng = 2*M_PI*kappaSq*qi;  // 2pi alpha^2*f_t for SS scheme exp(-r^2/alpha^2)

#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (j = 0; j < numparsS; j++) {
                tx = xi - xS[j];
                ty = yi - yS[j];
                tz = zi - zS[j];
                rad = sqrt(tx*tx + ty*ty + tz*tz);
                if (rad > 1e-14) {
                    teng = teng + ( qS[j] - qi * exp(-rad*rad/kappaSq)) * wS[j] / rad;
                }
            }
            denergy[i] += teng;
        }
#ifdef OPENACC_ENABLED
        }
#endif

    } else if (pot_type == 3 || pot_type == 7) {
#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (i = 0; i < numparsT; i++) {
            xi = xT[i];
            yi = yT[i];
            zi = zT[i];
            qi = qT[i];
            teng = 4*M_PI*qi/kappa/kappa;  // 4pi*f_t/k^2
            
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (j = 0; j < numparsS; j++) {
                tx = xi - xS[j];
                ty = yi - yS[j];
                tz = zi - zS[j];
                rad = sqrt(tx*tx + ty*ty + tz*tz);
                if (rad > 1e-14) {
                    teng += ( qS[j] - qi) * wS[j] * exp(-kappa * rad) / rad;
                }
            }
            denergy[i] += teng;
        }
        
#ifdef OPENACC_ENABLED
        }
#endif
    } // end pot=3 or 7
    
#ifdef OPENACC_ENABLED
    }
#endif

    *dpeng = sum(denergy, numparsT);

    return;
}
