#include <math.h>
#include <float.h>
#include <stdio.h>

#include "../../run_params/struct_run_params.h"
#include "sin-over-r_pc.h"


void K_SinOverR_PC_Lagrange(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
         int starting_index_of_target, int starting_index_of_cluster,
         double *target_x, double *target_y, double *target_z,
         double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
         struct RunParams *run_params, double *potential, int gpu_async_stream_id)
{

    double kernel_parameter = run_params->kernel_params[0];

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, \
                        cluster_x, cluster_y, cluster_z, cluster_charge, potential)
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
#ifdef OPENACC_ENABLED
            #pragma acc cache(cluster_x[starting_index_of_cluster : starting_index_of_cluster+number_of_interpolation_points_in_cluster], \
                              cluster_y[starting_index_of_cluster : starting_index_of_cluster+number_of_interpolation_points_in_cluster], \
                              cluster_z[starting_index_of_cluster : starting_index_of_cluster+number_of_interpolation_points_in_cluster], \
                              cluster_charge[starting_index_of_cluster : starting_index_of_cluster+number_of_interpolation_points_in_cluster])
#endif

            int jj = starting_index_of_cluster + j;
            double dx = tx - cluster_x[jj];
            double dy = ty - cluster_y[jj];
            double dz = tz - cluster_z[jj];
            double r  = sqrt(dx*dx + dy*dy + dz*dz);

            temporary_potential += cluster_charge[starting_index_of_cluster + j] * sin(kernel_parameter * r) / r;

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




void K_SinOverR_PC_Hermite(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster, int total_number_interpolation_points,
        double *target_x, double *target_y, double *target_z,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
        struct RunParams *run_params, double *potential, int gpu_async_stream_id)
{

    // total_number_interpolation_points is the stride, separating clustersQ, clustersQx, clustersQy, etc.
    double *cluster_charge_          = &cluster_charge[8*starting_index_of_cluster + 0*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_x   = &cluster_charge[8*starting_index_of_cluster + 1*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_y   = &cluster_charge[8*starting_index_of_cluster + 2*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_z   = &cluster_charge[8*starting_index_of_cluster + 3*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_xy  = &cluster_charge[8*starting_index_of_cluster + 4*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_yz  = &cluster_charge[8*starting_index_of_cluster + 5*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_xz  = &cluster_charge[8*starting_index_of_cluster + 6*number_of_interpolation_points_in_cluster];
    double *cluster_charge_delta_xyz = &cluster_charge[8*starting_index_of_cluster + 7*number_of_interpolation_points_in_cluster];

    double k  = run_params->kernel_params[0];
    double k2 = k * k;
    double k3 = k * k2;


#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x, target_y, target_z, \
                        cluster_x, cluster_y, cluster_z, cluster_charge, potential, \
                        cluster_charge_, cluster_charge_delta_x, cluster_charge_delta_y, cluster_charge_delta_z, \
                        cluster_charge_delta_xy, cluster_charge_delta_yz, cluster_charge_delta_xz, \
                        cluster_charge_delta_xyz)
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

                double rinv  = 1 / r;
                double r2inv = rinv  *  rinv;
                double r3inv = rinv  * r2inv;
                double r4inv = r2inv * r2inv;
                double r5inv = r3inv * r2inv;
                double r6inv = r3inv * r3inv;
                double r7inv = r4inv * r3inv;

                double sinr = sin(k*r);
                double cosr = cos(k*r);

                temporary_potential +=  sinr * rinv * (cluster_charge_[j])
                                    +  (sinr * r3inv  -  k * cosr * r2inv)
                                             * (cluster_charge_delta_x[j]*dx + cluster_charge_delta_y[j]*dy
                                              + cluster_charge_delta_z[j]*dz)
                                    +  (sinr * (3 * r5inv - k2 * r3inv)  -  3 * k * cosr * r4inv) 
                                             * (cluster_charge_delta_xy[j]*dx*dy + cluster_charge_delta_yz[j]*dy*dz
                                              + cluster_charge_delta_xz[j]*dx*dz)
                                    +  (sinr * (15 * r7inv - 6 * k2 * r5inv)  +  cosr * (k3 * r4inv - 15 * k * r6inv)) 
                                             *  cluster_charge_delta_xyz[j]*dx*dy*dz;

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
