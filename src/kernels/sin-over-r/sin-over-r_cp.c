#include <math.h>
#include <float.h>
#include <stdio.h>

#include "../../run_params/struct_run_params.h"
#include "sin-over-r_cp.h"


void K_SinOverR_CP_Lagrange(int number_of_sources_in_batch, int number_of_interpolation_points_in_cluster,
         int starting_index_of_sources, int starting_index_of_cluster,
         double *source_x, double *source_y, double *source_z, double *source_q,
         double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
         struct RunParams *run_params, int gpu_async_stream_id)
{

    double kernel_parameter = run_params->kernel_params[0];

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(source_x, source_y, source_z, source_q, \
                        cluster_x, cluster_y, cluster_z, cluster_q)
    {
#endif
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif	
    for (int i = 0; i < number_of_interpolation_points_in_cluster; i++) {

        double temporary_potential = 0.0;

        double cx = cluster_x[starting_index_of_cluster + i];
        double cy = cluster_y[starting_index_of_cluster + i];
        double cz = cluster_z[starting_index_of_cluster + i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential)
#endif
        for (int j = 0; j < number_of_sources_in_batch; j++) {
#ifdef OPENACC_ENABLED
            #pragma acc cache(source_x[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_y[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_z[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_q[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch])
#endif

            int jj = starting_index_of_sources + j;
            double dx = cx - source_x[jj];
            double dy = cy - source_y[jj];
            double dz = cz - source_z[jj];
            double r = sqrt(dx*dx + dy*dy + dz*dz);

            temporary_potential += source_q[jj] * sin(kernel_parameter * r) / r;

        } // end loop over interpolation points
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        cluster_q[starting_index_of_cluster + i] += temporary_potential;
    }
#ifdef OPENACC_ENABLED
    } // end kernel
#endif
    return;
}




void K_SinOverR_CP_Hermite(int number_of_sources_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_sources, int starting_index_of_cluster,
        double *source_x, double *source_y, double *source_z, double *source_q,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
        struct RunParams *run_params, int gpu_async_stream_id)
{

    double *cluster_q_     = &cluster_q[8*starting_index_of_cluster + 0*number_of_interpolation_points_in_cluster];
    double *cluster_q_dx   = &cluster_q[8*starting_index_of_cluster + 1*number_of_interpolation_points_in_cluster];
    double *cluster_q_dy   = &cluster_q[8*starting_index_of_cluster + 2*number_of_interpolation_points_in_cluster];
    double *cluster_q_dz   = &cluster_q[8*starting_index_of_cluster + 3*number_of_interpolation_points_in_cluster];
    double *cluster_q_dxy  = &cluster_q[8*starting_index_of_cluster + 4*number_of_interpolation_points_in_cluster];
    double *cluster_q_dyz  = &cluster_q[8*starting_index_of_cluster + 5*number_of_interpolation_points_in_cluster];
    double *cluster_q_dxz  = &cluster_q[8*starting_index_of_cluster + 6*number_of_interpolation_points_in_cluster];
    double *cluster_q_dxyz = &cluster_q[8*starting_index_of_cluster + 7*number_of_interpolation_points_in_cluster];

    double k  = run_params->kernel_params[0];
    double k2 = k * k;
    double k3 = k * k2;

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(source_x, source_y, source_z, source_q, \
                        cluster_x, cluster_y, cluster_z, \
                        cluster_q_, cluster_q_dx, cluster_q_dy, cluster_q_dz, \
                        cluster_q_dxy, cluster_q_dyz, cluster_q_dxz, \
                        cluster_q_dxyz)
    {
#endif
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < number_of_interpolation_points_in_cluster; i++) {

        double temp_pot_     = 0.0;
        double temp_pot_dx   = 0.0;
        double temp_pot_dy   = 0.0;
        double temp_pot_dz   = 0.0;
        double temp_pot_dxy  = 0.0;
        double temp_pot_dyz  = 0.0;
        double temp_pot_dxz  = 0.0;
        double temp_pot_dxyz = 0.0;
        
        int ii = starting_index_of_cluster + i;
        double cx = cluster_x[ii];
        double cy = cluster_y[ii];
        double cz = cluster_z[ii];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp_pot_dx)  reduction(+:temp_pot_dy)  reduction(+:temp_pot_dz)  \
                                     reduction(+:temp_pot_dxy) reduction(+:temp_pot_dyz) reduction(+:temp_pot_dxz) \
                                     reduction(+:temp_pot_)    reduction(+:temp_pot_dxyz)
                                               
#endif
        for (int j = 0; j < number_of_sources_in_batch; j++) {
#ifdef OPENACC_ENABLED
            #pragma acc cache(source_x[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_y[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_z[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_q[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch])
#endif

            int jj = starting_index_of_sources + j;
            double dx = source_x[jj] - cx;
            double dy = source_y[jj] - cy;
            double dz = source_z[jj] - cz;
            double r = sqrt(dx*dx + dy*dy + dz*dz);
            

            double rinv  = 1 / r;
            double r2inv = rinv  *  rinv;
            double r3inv = rinv  * r2inv;
            double r4inv = r2inv * r2inv;
            double r5inv = r3inv * r2inv;
            double r6inv = r3inv * r3inv;
            double r7inv = r4inv * r3inv;

            double sinr = sin(k*r) * source_q[jj];
            double cosr = cos(k*r) * source_q[jj];

            double term_d0 = sinr * rinv;
            double term_d1 = sinr * r3inv  -  k * cosr * r2inv;
            double term_d2 = sinr * (3 * r5inv - k2 * r3inv)  -  3 * k * cosr * r4inv;
            double term_d3 = sinr * (15 * r7inv - 6 * k2 * r5inv)  +  cosr * (k3 * r4inv - 15 * k * r6inv);

            temp_pot_     += term_d0;
            temp_pot_dx   += term_d1 * dx;
            temp_pot_dy   += term_d1 * dy;
            temp_pot_dz   += term_d1 * dz;
            temp_pot_dxy  += term_d2 * dx * dy;
            temp_pot_dyz  += term_d2 * dy * dz;
            temp_pot_dxz  += term_d2 * dx * dz;
            temp_pot_dxyz += term_d3 * dx * dy * dz;


        } // end loop over interpolation points
        
#ifdef OPENACC_ENABLED
        #pragma acc atomic
        cluster_q_[i]     += temp_pot_;
        #pragma acc atomic
        cluster_q_dx[i]   += temp_pot_dx;
        #pragma acc atomic
        cluster_q_dy[i]   += temp_pot_dy;
        #pragma acc atomic
        cluster_q_dz[i]   += temp_pot_dz;
        #pragma acc atomic
        cluster_q_dxy[i]  += temp_pot_dxy;
        #pragma acc atomic
        cluster_q_dyz[i]  += temp_pot_dyz;
        #pragma acc atomic
        cluster_q_dxz[i]  += temp_pot_dxz;
        #pragma acc atomic
        cluster_q_dxyz[i] += temp_pot_dxyz;
#else
        cluster_q_[i]     += temp_pot_;
        cluster_q_dx[i]   += temp_pot_dx;
        cluster_q_dy[i]   += temp_pot_dy;
        cluster_q_dz[i]   += temp_pot_dz;
        cluster_q_dxy[i]  += temp_pot_dxy;
        cluster_q_dyz[i]  += temp_pot_dyz;
        cluster_q_dxz[i]  += temp_pot_dxz;
        cluster_q_dxyz[i] += temp_pot_dxyz;
#endif


    }
#ifdef OPENACC_ENABLED
    } // end kernel
#endif
    return;
}
