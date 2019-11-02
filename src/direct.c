#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <float.h>

#include "array.h"
#include "tools.h"
#include "direct.h"



void direct_eng(double *xS, double *yS, double *zS, double *qS, double *wS,
				double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng, double kappa,
                double (*kernel)(double,  double,  double,  double,  double,  double,  double,  double,  double, double))
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
            teng = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (j = 0; j < numparsS; j++) {
//                tx = xi - xS[j];
//                ty = yi - yS[j];
//                tz = zi - zS[j];
//                rad = sqrt(tx*tx + ty*ty + tz*tz);
//                if (rad > 1e-14) {
//                    teng = teng + qS[j] * wS[j] / rad;
				teng = teng + kernel(xi, yi, zi, qi, xS[j], yS[j], zS[j], qS[j], wS[j], kappa);
//                }
            }
            denergy[i] +=  teng;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif
#ifdef OPENACC_ENABLED
    } // end acc data region
#endif

    *dpeng = sum(denergy, numparsT);

    return;
}
