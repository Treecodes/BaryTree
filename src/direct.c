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



void directSummation(double *source_x, double *source_y, double *source_z, double *source_charge, double *source_weight,
                double *target_x, double *target_y, double *target_z, double *target_charge,
                int number_of_sources, int number_of_targets, double *potential, double *total_potential, double kernel_parameter,
                char *kernel_name)
{

#ifdef OPENACC_ENABLED
    #pragma acc data copyin (source_x[0:number_of_sources], source_y[0:number_of_sources], source_z[0:number_of_sources], \
                             source_charge[0:number_of_sources], source_weight[0:number_of_sources], source_x[0:number_of_targets], \
                             source_y[0:number_of_targets], source_z[0:number_of_targets], source_charge[0:number_of_targets])
    {
#endif



/***************************************/
/************* Coulomb *****************/
/***************************************/

    if (strcmp(kernel_name, "coulomb") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < number_of_targets; i++) {

            double temporary_potential = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < number_of_sources; j++)
                temporary_potential += coulombKernel(target_x[i], target_y[i], target_z[i], target_charge[i], source_x[j], source_y[j], source_z[j], source_charge[j], source_weight[j], kernel_parameter);

            potential[i] += temporary_potential;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif



/***************************************/
/*************** Yukawa ****************/
/***************************************/

    } else if (strcmp(kernel_name, "yukawa") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < number_of_targets; i++) {

            double temporary_potential = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < number_of_sources; j++)
                temporary_potential += yukawaKernel(target_x[i], target_y[i], target_z[i], target_charge[i], source_x[j], source_y[j], source_z[j], source_charge[j], source_weight[j], kernel_parameter);

            potential[i] += temporary_potential;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif



/***************************************/
/********** Coulomb with SS ************/
/***************************************/

    } else if (strcmp(kernel_name, "coulomb_SS") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < number_of_targets; i++) {

            double temporary_potential = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < number_of_sources; j++)
                temporary_potential += coulombKernel_SS_direct(target_x[i], target_y[i], target_z[i], target_charge[i], source_x[j], source_y[j], source_z[j], source_charge[j], source_weight[j], kernel_parameter);

            potential[i] += temporary_potential;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif



/***************************************/
/********** Yukawa with SS *************/
/***************************************/

    } else if (strcmp(kernel_name, "yukawa_SS") == 0) {

#ifdef OPENACC_ENABLED
        #pragma acc kernels
        {
        #pragma acc loop independent
#endif
        for (int i = 0; i < number_of_targets; i++) {

            double temporary_potential = 0.0;
                
#ifdef OPENACC_ENABLED
            #pragma acc loop independent
#endif
            for (int j = 0; j < number_of_sources; j++)
                temporary_potential += yukawaKernel_SS_direct(target_x[i], target_y[i], target_z[i], target_charge[i], source_x[j], source_y[j], source_z[j], source_charge[j], source_weight[j], kernel_parameter);

            potential[i] += temporary_potential;
        }
#ifdef OPENACC_ENABLED
        } // end acc kernels
#endif

    } // end kernel selection


#ifdef OPENACC_ENABLED
    } // end acc data region
#endif

    *total_potential = sum(potential, number_of_targets);

    return;
}
