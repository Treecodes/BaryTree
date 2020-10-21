#ifdef OPENACC_ENABLED
    #include <accelmath.h>
    #define M_PI 3.14159265358979323846264338327950288
#else
    #include <math.h>
#endif
#include <stdio.h>

#include "../../run_params/struct_run_params.h"
#include "mq_pp.h"


void K_MQ_PP(int number_of_targets_in_batch, int number_of_source_points_in_cluster,
        int starting_index_of_target, int starting_index_of_source,
        double *target_x, double *target_y, double *target_z,
        double *source_x, double *source_y, double *source_z, double *source_charge,
        struct RunParams *run_params, double *potential, int gpu_async_stream_id)
{

    double domainLength = run_params->kernel_params[0];
    double delta = run_params->kernel_params[1];
    double deltaLsq = delta * delta / domainLength / domainLength;
    double norm_delta_L = sqrt(1 + 4 * deltaLsq);
    

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, \
                        source_x, source_y, source_z, source_charge, potential)
    {
#endif
#ifdef OPENACC_ENABLED
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
            double dz = (tz - source_z[jj]) / domainLength;

            if (dz < -0.5) {
                dz += 1.0;
            }
            if (dz > 0.5) {
                dz -= 1.0;
            }
            temporary_potential += source_charge[jj]
                                * (.5 * dz * norm_delta_L / sqrt(dz * dz + deltaLsq) - dz);
        } // end loop over interpolation points
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        potential[ii] += temporary_potential;
    }
#ifdef OPENACC_ENABLED
    } // end kernel
#endif
    return;
}
