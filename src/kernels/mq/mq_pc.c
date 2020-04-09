#ifdef OPENACC_ENABLED
    #include <accelmath.h>
    #define M_PI 3.14159265358979323846264338327950288
#else
    #include <math.h>
#endif
#include <stdio.h>

#include "../../run_params/struct_run_params.h"
#include "mq_pc.h"


void K_MQ_PC_Lagrange(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster,
        double *target_x, double *target_y, double *target_z,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
        struct RunParams *run_params, double *potential, int gpu_async_stream_id)
{

    double domainLength = run_params->kernel_params[0];
    double delta = run_params->kernel_params[1];
    double deltaLsq = delta * delta / domainLength / domainLength;
    double norm_delta_L = sqrt(1 + 4 * deltaLsq);

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, \
                        cluster_x, cluster_y, cluster_z, cluster_charge, potential)
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
        for (int j = 0; j < number_of_interpolation_points_in_cluster; j++) {

            int jj = starting_index_of_cluster + j;
            double dz = (tz - cluster_z[jj]) / domainLength;

            if (dz < -0.5) {
                dz += 1.0;
            }
            if (dz > 0.5) {
                dz -= 1.0;
            }

            temporary_potential += cluster_charge[jj]
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




void K_MQ_PC_Hermite(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster, int total_number_interpolation_points,
        double *target_x, double *target_y, double *target_z,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
        struct RunParams *run_params, double *potential, int gpu_async_stream_id)
{
    printf("[BaryTree] ERROR! MQ KERNEL NOT IMPLEMENTED FOR HERMITE. Exiting.\n");
    return;
}
