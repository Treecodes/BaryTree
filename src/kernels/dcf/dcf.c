#include <math.h>
#include <float.h>
#include <stdio.h>

#include "dcf.h"

void dcfDirect(int number_of_targets_in_batch, int number_of_source_points_in_cluster,
        int starting_index_of_target, int starting_index_of_source,
        double *target_x, double *target_y, double *target_z,
        double *source_x, double *source_y, double *source_z, double *source_q, double *source_w,
        double kernel_parameter, double *potential, int gpu_async_stream_id)
{

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
                temporary_potential += source_q[starting_index_of_source + j]
                                     * source_w[starting_index_of_source + j]
                                     * erf(r / kernel_parameter) / r;
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




void dcfApproximationLagrange(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster,
        double *target_x, double *target_y, double *target_z,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
        double kernel_parameter, double *potential, int gpu_async_stream_id)
{

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

            if (r > DBL_MIN){
                temporary_potential += cluster_q[starting_index_of_cluster + j]
                                     * erf(r / kernel_parameter) / r;
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




void dcfApproximationHermite(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster, int total_number_interpolation_points,
        double *target_x, double *target_y, double *target_z,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
        double kernel_parameter, double *potential, int gpu_async_stream_id)
{

    // total_number_interpolation_points is the stride, separating clustersQ, clustersQx, clustersQy, etc.
    double *cluster_q_     = &cluster_q[8*starting_index_of_cluster + 0*number_of_interpolation_points_in_cluster];
    double *cluster_q_dx   = &cluster_q[8*starting_index_of_cluster + 1*number_of_interpolation_points_in_cluster];
    double *cluster_q_dy   = &cluster_q[8*starting_index_of_cluster + 2*number_of_interpolation_points_in_cluster];
    double *cluster_q_dz   = &cluster_q[8*starting_index_of_cluster + 3*number_of_interpolation_points_in_cluster];
    double *cluster_q_dxy  = &cluster_q[8*starting_index_of_cluster + 4*number_of_interpolation_points_in_cluster];
    double *cluster_q_dyz  = &cluster_q[8*starting_index_of_cluster + 5*number_of_interpolation_points_in_cluster];
    double *cluster_q_dxz  = &cluster_q[8*starting_index_of_cluster + 6*number_of_interpolation_points_in_cluster];
    double *cluster_q_dxyz = &cluster_q[8*starting_index_of_cluster + 7*number_of_interpolation_points_in_cluster];

    double twoinvEta2 = 2.0 / (kernel_parameter * kernel_parameter);
    double teninvEta2 = 5.0 * twoinvEta2;
    double fourinvEta4 = twoinvEta2 * twoinvEta2;

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, \
                        cluster_x, cluster_y, cluster_z, potential, \
                        cluster_q_, cluster_q_dx, cluster_q_dy, cluster_q_dz, \
                        cluster_q_dxy, cluster_q_dyz, cluster_q_dxz, \
                        cluster_q_dxyz)
    {
#endif
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < number_of_targets_in_batch; i++) {

        int ii = starting_index_of_target + i;
        double temporary_potential = 0.0;

        double tx = target_x[ii];
        double ty = target_y[ii];
        double tz = target_z[ii];

        #ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential)
#endif
        for (int j = 0; j < number_of_interpolation_points_in_cluster; j++) {

            int jj = starting_index_of_cluster + j;
            double dx = tx - cluster_x[jj];
            double dy = ty - cluster_y[jj];
            double dz = tz - cluster_z[jj];
            double r  = sqrt(dx*dx + dy*dy + dz*dz);

            if (r > DBL_MIN) {

                double rinv = 1 / r;
                double r2inv = rinv * rinv;
                double r4inv = r2inv * r2inv;
                double r_eta = r / kernel_parameter;

                double pot = erf(r_eta) * rinv; 
                double auxpot = 2.0 / sqrt(M_PI) * exp(-r_eta * r_eta);
                double auxpot_eta = auxpot / kernel_parameter;

                double dpot1 = r2inv * (pot - auxpot);
                double dpot2 = r2inv * ( 3.0 * pot * r2inv - auxpot_eta
                             * ( 3.0 * r2inv + twoinvEta2));
                double dpot3 = r2inv * (15.0 * pot * r4inv - auxpot_eta
                             * (15.0 * r4inv + teninvEta2 * r2inv + fourinvEta4));

                temporary_potential +=  pot  *     cluster_q_[j] 
                                     + dpot1 *  (cluster_q_dx[j] * dx
                                             +   cluster_q_dy[j] * dy
                                             +   cluster_q_dz[j] * dz)
                                     + dpot2 * (cluster_q_dxy[j] * dx * dy
                                             +  cluster_q_dyz[j] * dy * dz
                                             +  cluster_q_dxz[j] * dx * dz)
                                     + dpot3 * cluster_q_dxyz[j] * dx * dy * dz;

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
