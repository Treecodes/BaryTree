#include <math.h>
#include <float.h>
#include <stdio.h>

#include "../struct_kernel.h"

#include "yukawa_singularity_subtraction.h"

void yukawaSingularitySubtractionDirect(int number_of_targets_in_batch, int number_of_source_points_in_cluster,
        int starting_index_of_target, int starting_index_of_source,
        double *target_x, double *target_y, double *target_z, double *target_charge,
        double *source_x, double *source_y, double *source_z, double *source_charge, double * source_weight,
        struct kernel *kernel, double *potential, int gpu_async_stream_id)
{

    double kernel_parameter=kernel->parameters[0];

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, target_charge, \
                        source_x, source_y, source_z, source_charge, source_weight, potential)
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
        double tq = target_charge[ii];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential)
#endif
        for (int j = 0; j < number_of_source_points_in_cluster; j++) {

            int jj = starting_index_of_source + j;
            double dx = tx - source_x[jj];
            double dy = ty - source_y[jj];
            double dz = tz - source_z[jj];
            double r  = sqrt(dx*dx + dy*dy + dz*dz);

            if (r > DBL_MIN) {
                temporary_potential += (source_charge[jj] - tq) * source_weight[jj] * exp(-kernel_parameter*r) / r;
            }
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


void yukawaSingularitySubtractionApproximationLagrange(int number_of_targets_in_batch,
        int number_of_interpolation_points_in_cluster, int starting_index_of_target, int starting_index_of_cluster,
        double *target_x, double *target_y, double *target_z, double *target_charge,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge, double *cluster_weight,
        struct kernel *kernel, double *potential, int gpu_async_stream_id)
{
    double kernel_parameter=kernel->parameters[0];

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, target_charge, \
                        cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight, potential)
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
        double tq = target_charge[ii];

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
                temporary_potential += (cluster_charge[jj] - tq * cluster_weight[jj]) * exp(-kernel_parameter*r) /r;
            }
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


void yukawaSingularitySubtractionApproximationHermite(int number_of_targets_in_batch,
        int number_of_interpolation_points_in_cluster, int starting_index_of_target,
        int starting_index_of_cluster, int total_number_interpolation_points,
        double *target_x, double *target_y, double *target_z, double *target_charge,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge, double *cluster_weight,
        struct kernel *kernel, double *potential, int gpu_async_stream_id)
{

    double kernel_parameter=kernel->parameters[0];
    double kernel_parameter2 = kernel_parameter * kernel_parameter;
    double kernel_parameter3 = kernel_parameter * kernel_parameter2;

    // total_number_interpolation_points is the stride, separating clustersQ, clustersQx, clustersQy, etc.
    double *cluster_charge_          = &cluster_charge[8*starting_index_of_cluster + 0*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_x   = &cluster_charge[8*starting_index_of_cluster + 1*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_y   = &cluster_charge[8*starting_index_of_cluster + 2*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_z   = &cluster_charge[8*starting_index_of_cluster + 3*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_xy  = &cluster_charge[8*starting_index_of_cluster + 4*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_yz  = &cluster_charge[8*starting_index_of_cluster + 5*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_xz  = &cluster_charge[8*starting_index_of_cluster + 6*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_xyz = &cluster_charge[8*starting_index_of_cluster + 7*number_of_interpolation_points_in_cluster];

    double *cluster_weight_          = &cluster_weight[8*starting_index_of_cluster + 0*number_of_interpolation_points_in_cluster];
    double *cluster_weight_delta_x   = &cluster_weight[8*starting_index_of_cluster + 1*number_of_interpolation_points_in_cluster];
    double *cluster_weight_delta_y   = &cluster_weight[8*starting_index_of_cluster + 2*number_of_interpolation_points_in_cluster];
    double *cluster_weight_delta_z   = &cluster_weight[8*starting_index_of_cluster + 3*number_of_interpolation_points_in_cluster];
    double *cluster_weight_delta_xy  = &cluster_weight[8*starting_index_of_cluster + 4*number_of_interpolation_points_in_cluster];
    double *cluster_weight_delta_yz  = &cluster_weight[8*starting_index_of_cluster + 5*number_of_interpolation_points_in_cluster];
    double *cluster_weight_delta_xz  = &cluster_weight[8*starting_index_of_cluster + 6*number_of_interpolation_points_in_cluster];
    double *cluster_weight_delta_xyz = &cluster_weight[8*starting_index_of_cluster + 7*number_of_interpolation_points_in_cluster];


#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, target_charge, \
                        cluster_x, cluster_y, cluster_z, cluster_charge, cluster_weight, potential, \
                        cluster_charge_, cluster_charge_delta_x, cluster_charge_delta_y, cluster_charge_delta_z, \
                        cluster_charge_delta_xy, cluster_charge_delta_yz, cluster_charge_delta_xz, \
                        cluster_charge_delta_xyz, \
                        cluster_weight_, cluster_weight_delta_x, cluster_weight_delta_y, cluster_weight_delta_z, \
                        cluster_weight_delta_xy, cluster_weight_delta_yz, cluster_weight_delta_xz, \
                        cluster_weight_delta_xyz)
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
        double tq = target_charge[ii];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential)
#endif
        for (int j = 0; j < number_of_interpolation_points_in_cluster; j++) {

            int jj = starting_index_of_cluster + j;
            double dx = tx - cluster_x[jj];
            double dy = ty - cluster_y[jj];
            double dz = tz - cluster_z[jj];
            double r  = sqrt( dx*dx + dy*dy + dz*dz);

            double r2 = r*r;
            double r3 = r2*r;
            double r4 = r2*r2;
            double rinv = 1 / r;
            double r3inv = rinv*rinv*rinv;
            double r5inv = r3inv*rinv*rinv;
            double r7inv = r5inv*rinv*rinv;

            double kr = kernel_parameter * r;
            double k2r2 = kr * kr;
            double k3r3 = k2r2 * kr;

            double charge_diff    = cluster_charge_[j]          - cluster_weight_[j]          * tq;
            double delta_x_diff   = cluster_charge_delta_x[j]   - cluster_weight_delta_x[j]   * tq;
            double delta_y_diff   = cluster_charge_delta_y[j]   - cluster_weight_delta_y[j]   * tq;
            double delta_z_diff   = cluster_charge_delta_z[j]   - cluster_weight_delta_z[j]   * tq;
            double delta_xy_diff  = cluster_charge_delta_xy[j]  - cluster_weight_delta_xy[j]  * tq;
            double delta_yz_diff  = cluster_charge_delta_yz[j]  - cluster_weight_delta_yz[j]  * tq;
            double delta_xz_diff  = cluster_charge_delta_xz[j]  - cluster_weight_delta_xz[j]  * tq;
            double delta_xyz_diff = cluster_charge_delta_xyz[j] - cluster_weight_delta_xyz[j] * tq;

            if (r > DBL_MIN) {

                temporary_potential += exp(-kernel_parameter*r)
                              * (rinv  * (charge_diff)
                               + r3inv * (1 + kernel_parameter*r)
                                       * (delta_x_diff*dx + delta_y_diff*dy + delta_z_diff*dz)
                               + r5inv * (3 + 3*kernel_parameter*r + kernel_parameter2*r2)
                                       * (delta_xy_diff*dx*dy + delta_yz_diff*dy*dz + delta_xz_diff*dx*dz)
                               + r7inv * (15 + 15*kernel_parameter*r + 6*kernel_parameter2*r2 + kernel_parameter3*r3)
                                       * delta_xyz_diff*dx*dy*dz);

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



void yukawaSingularitySubtractionCorrection(double *potential, double *target_q,
                                            int numTargets, struct kernel *kernel)
{
    double kernel_parameter=kernel->parameters[0];
    double param = 4.0 * M_PI / kernel_parameter / kernel_parameter;
    for (int i = 0; i < numTargets; i++) potential[i] += param * target_q[i];

    return;
}
