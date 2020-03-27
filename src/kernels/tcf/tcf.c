#include <math.h>
#include <float.h>
#include <stdio.h>

#include "tcf.h"

void tcfDirect(int number_of_targets_in_batch, int number_of_source_points_in_cluster,
        int starting_index_of_target, int starting_index_of_source,
        double *target_x, double *target_y, double *target_z,
        double *source_x, double *source_y, double *source_z, double *source_q, double *source_w,
        double kernel_parameter1, double kernel_parameter2, double *potential, int gpu_async_stream_id)
{

    double kap_eta_2 = kernel_parameter1 * kernel_parameter2 / 2.0;

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, \
                        source_x, source_y, source_z, source_q, source_w, potential)
    {
#endif
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
        for (int i = 0; i < number_of_targets_in_batch; i++) {

        double temporary_potential = 0.0;

        double tx = target_x[starting_index_of_target + i];
        double ty = target_y[starting_index_of_target + i];
        double tz = target_z[starting_index_of_target + i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential)
#endif
        for (int j = 0; j < number_of_source_points_in_cluster; j++) {

            double dx = tx - source_x[starting_index_of_source + j];
            double dy = ty - source_y[starting_index_of_source + j];
            double dz = tz - source_z[starting_index_of_source + j];
            double r  = sqrt(dx*dx + dy*dy + dz*dz);

            if (r > DBL_MIN) {
                double kap_r = kernel_parameter1 * r;
                double r_eta = r / kernel_parameter2;
                temporary_potential += source_q[starting_index_of_source + j]
                                     * source_w[starting_index_of_source + j] / r
                                     * (exp(-kap_r) * erfc(kap_eta_2 - r_eta)
                                     -  exp( kap_r) * erfc(kap_eta_2 + r_eta));
            }
        } // end loop over interpolation points
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        potential[starting_index_of_target + i] += temporary_potential;
    }
#ifdef OPENACC_ENABLED
    } // end kernel
#endif
    return;
}




void tcfApproximationLagrange(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster,
        double *target_x, double *target_y, double *target_z,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
        double kernel_parameter1, double kernel_parameter2, double *potential, int gpu_async_stream_id)
{

    double kap_eta_2 = kernel_parameter1 * kernel_parameter2 / 2.0;

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, \
                        cluster_x, cluster_y, cluster_z, cluster_q, potential)
    {
#endif
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < number_of_targets_in_batch; i++) {

        double temporary_potential = 0.0;

        double tx = target_x[starting_index_of_target + i];
        double ty = target_y[starting_index_of_target + i];
        double tz = target_z[starting_index_of_target + i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential)
#endif
        for (int j = 0; j < number_of_interpolation_points_in_cluster; j++) {

            double dx = tx - cluster_x[starting_index_of_cluster + j];
            double dy = ty - cluster_y[starting_index_of_cluster + j];
            double dz = tz - cluster_z[starting_index_of_cluster + j];
            double r  = sqrt(dx*dx + dy*dy + dz*dz);

            if (r > DBL_MIN) {
                double kap_r = kernel_parameter1 * r;
                double r_eta = r / kernel_parameter2;
                temporary_potential += cluster_q[starting_index_of_cluster + j] / r
                                     * (exp(-kap_r) * erfc(kap_eta_2 - r_eta)
                                     -  exp( kap_r) * erfc(kap_eta_2 + r_eta));
            }
        } // end loop over interpolation points
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        potential[starting_index_of_target + i] += temporary_potential;
    }
#ifdef OPENACC_ENABLED
    } // end kernel
#endif
    return;
}
