#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <float.h>

#include "array.h"
#include "tools.h"

#include "kernels/kernels.h"
#include "direct.h"



void direct_eng(double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng, double kappa,
                char *kernelName)
{

#ifdef OPENACC_ENABLED
    #pragma acc data copyin (xS[0:numparsS], yS[0:numparsS], zS[0:numparsS], \
                             qS[0:numparsS], wS[0:numparsS], xT[0:numparsT], \
                             yT[0:numparsT], zT[0:numparsT], qT[0:numparsT])
    {
#endif



/***************************************/
/************* Coulomb *****************/
/***************************************/

    if (strcmp(kernelName, "coulomb") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < numparsT; i++) {
            double xi = xT[i];
            double yi = yT[i];
            double zi = zT[i];
            double qi = qT[i];
            double teng = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < numparsS; j++)
                teng += coulombKernel(xi, yi, zi, qi, xS[j], yS[j], zS[j], qS[j], wS[j], kappa);

            denergy[i] += teng;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif



/***************************************/
/*************** Yukawa ****************/
/***************************************/

    } else if (strcmp(kernelName, "yukawa") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < numparsT; i++) {
            double xi = xT[i];
            double yi = yT[i];
            double zi = zT[i];
            double qi = qT[i];
            double teng = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < numparsS; j++)
                teng += yukawaKernel(xi, yi, zi, qi, xS[j], yS[j], zS[j], qS[j], wS[j], kappa);

            denergy[i] += teng;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif



/***************************************/
/********** Coulomb with SS ************/
/***************************************/

    } else if (strcmp(kernelName, "coulomb_SS") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < numparsT; i++) {
            double xi = xT[i];
            double yi = yT[i];
            double zi = zT[i];
            double qi = qT[i];
            double teng = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < numparsS; j++)
                teng += coulombKernel_SS_direct(xi, yi, zi, qi, xS[j], yS[j], zS[j], qS[j], wS[j], kappa);

            denergy[i] += teng;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif



/***************************************/
/********** Yukawa with SS *************/
/***************************************/

    } else if (strcmp(kernelName, "yukawa_SS") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < numparsT; i++) {
            double xi = xT[i];
            double yi = yT[i];
            double zi = zT[i];
            double qi = qT[i];
            double teng = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < numparsS; j++)
                teng += yukawaKernel_SS_direct(xi, yi, zi, qi, xS[j], yS[j], zS[j], qS[j], wS[j], kappa);

            denergy[i] += teng;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif

    } // end kernel selection


#ifdef OPENACC_ENABLED
    } // end acc data region
#endif

    *dpeng = sum(denergy, numparsT);

    return;
}
