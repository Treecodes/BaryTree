#include <math.h>
#include <float.h>
#include <stdio.h>

#include "../../struct_run_params.h"
#include "atan_direct.h"


void K_Atan_Direct(int number_of_targets_in_batch, int number_of_source_points_in_cluster,
        int starting_index_of_target, int starting_index_of_source,
        double *target_x, double *target_y, double *target_z,
        double *source_x, double *source_y, double *source_z, double *source_charge, double *source_weight,
        struct RunParams *run_params, double *potential, int gpu_async_stream_id)
{

    double domainLength=run_params->kernel_params[0];
    double delta=run_params->kernel_params[1];
    double wadj = 1/(1-delta/sqrt(1+delta*delta));

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, \
                        source_x, source_y, source_z, source_charge, source_weight, potential)
    {
    #pragma acc loop independent
#endif
    for (int i = 0; i < number_of_targets_in_batch; i++) {

        int ii = starting_index_of_target + i;
        double temporary_potential = 0.0;
        double tz = target_z[ii];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential)
#endif
        for (int j = 0; j < number_of_source_points_in_cluster; j++) {


            int jj = starting_index_of_source + j;
            double dz = (tz - source_z[jj])/domainLength;


            if (fabs(dz-0.5) > DBL_MIN) {
                if (fabs(dz+0.5) > DBL_MIN) {
//                temporary_potential += source_charge[jj] * source_weight[jj] * ( 1/M_PI * atan( sqrt( 1 + 1.0/(delta*delta))* tan(M_PI * dz)) - fmod(dz-0.5,1.0) + 0.5);
                }
            }


        } // end loop over interpolation points
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        potential[ii] += wadj * temporary_potential;
    }
#ifdef OPENACC_ENABLED
    } // end kernel
#endif
    return;
}
