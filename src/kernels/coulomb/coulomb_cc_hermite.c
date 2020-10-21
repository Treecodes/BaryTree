#include <math.h>
#include <float.h>
#include <stdio.h>

#include "../../run_params/struct_run_params.h"
#include "coulomb_cc_hermite.h"


void K_Coulomb_CC_Hermite(int number_of_sources_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_source_cluster, int starting_index_of_target_cluster,
        double *source_cluster_x, double *source_cluster_y, double *source_cluster_z, double *source_cluster_q,
        double *source_cluster_w,
        double *target_cluster_x, double *target_cluster_y, double *target_cluster_z, double *target_cluster_q,
        struct RunParams *run_params, int gpu_async_stream_id)
{

    double *source_cluster_q_     = &source_cluster_q[8*starting_index_of_source_cluster + 0*number_of_interpolation_points_in_cluster];
    double *source_cluster_q_dx   = &source_cluster_q[8*starting_index_of_source_cluster + 1*number_of_interpolation_points_in_cluster];
    double *source_cluster_q_dy   = &source_cluster_q[8*starting_index_of_source_cluster + 2*number_of_interpolation_points_in_cluster];
    double *source_cluster_q_dz   = &source_cluster_q[8*starting_index_of_source_cluster + 3*number_of_interpolation_points_in_cluster];
    double *source_cluster_q_dxy  = &source_cluster_q[8*starting_index_of_source_cluster + 4*number_of_interpolation_points_in_cluster];
    double *source_cluster_q_dyz  = &source_cluster_q[8*starting_index_of_source_cluster + 5*number_of_interpolation_points_in_cluster];
    double *source_cluster_q_dxz  = &source_cluster_q[8*starting_index_of_source_cluster + 6*number_of_interpolation_points_in_cluster];
    double *source_cluster_q_dxyz = &source_cluster_q[8*starting_index_of_source_cluster + 7*number_of_interpolation_points_in_cluster];

    double *target_cluster_q_     = &target_cluster_q[8*starting_index_of_target_cluster + 0*number_of_interpolation_points_in_cluster];
    double *target_cluster_q_dx   = &target_cluster_q[8*starting_index_of_target_cluster + 1*number_of_interpolation_points_in_cluster];
    double *target_cluster_q_dy   = &target_cluster_q[8*starting_index_of_target_cluster + 2*number_of_interpolation_points_in_cluster];
    double *target_cluster_q_dz   = &target_cluster_q[8*starting_index_of_target_cluster + 3*number_of_interpolation_points_in_cluster];
    double *target_cluster_q_dxy  = &target_cluster_q[8*starting_index_of_target_cluster + 4*number_of_interpolation_points_in_cluster];
    double *target_cluster_q_dyz  = &target_cluster_q[8*starting_index_of_target_cluster + 5*number_of_interpolation_points_in_cluster];
    double *target_cluster_q_dxz  = &target_cluster_q[8*starting_index_of_target_cluster + 6*number_of_interpolation_points_in_cluster];
    double *target_cluster_q_dxyz = &target_cluster_q[8*starting_index_of_target_cluster + 7*number_of_interpolation_points_in_cluster];



#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(source_cluster_x, source_cluster_y, source_cluster_z, \
                        source_cluster_q_,     source_cluster_q_dx,  source_cluster_q_dy,  source_cluster_q_dz, \
                        source_cluster_q_dxy,  source_cluster_q_dyz, source_cluster_q_dxz, \
                        source_cluster_q_dxyz, source_cluster_w, \
                        target_cluster_x,      target_cluster_y,     target_cluster_z, \
                        target_cluster_q_,     target_cluster_q_dx,  target_cluster_q_dy,  target_cluster_q_dz, \
                        target_cluster_q_dxy,  target_cluster_q_dyz, target_cluster_q_dxz, \
                        target_cluster_q_dxyz)
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
        
        int ii = starting_index_of_target_cluster + i;
        double cx = target_cluster_x[ii];
        double cy = target_cluster_y[ii];
        double cz = target_cluster_z[ii];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temp_pot_dx)  reduction(+:temp_pot_dy)  reduction(+:temp_pot_dz)  \
                                     reduction(+:temp_pot_dxy) reduction(+:temp_pot_dyz) reduction(+:temp_pot_dxz) \
                                     reduction(+:temp_pot_)    reduction(+:temp_pot_dxyz)
                                               
#endif
        for (int j = 0; j < number_of_sources_in_batch; j++) {
//#ifdef OPENACC_ENABLED
//            #pragma acc cache(source_x[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
//                              source_y[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
//                              source_z[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
//                              source_q[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
//                              source_w[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch])
//#endif

            int jj = starting_index_of_source_cluster + j;
            double dx = source_cluster_x[jj] - cx;
            double dy = source_cluster_y[jj] - cy;
            double dz = source_cluster_z[jj] - cz;
            double r2 = dx*dx + dy*dy + dz*dz;

            double r2inv = 1 / r2;
            double rinv  = 1 / sqrt(r2);
            double r3inv = rinv  * r2inv;
            double r5inv = r3inv * r2inv;
            double r7inv = r5inv * r2inv;
            
            r5inv *= 3.0;

            temp_pot_     += source_cluster_q_[j]     * rinv;
            temp_pot_dx   += source_cluster_q_dx[j]   * r3inv * dx;
            temp_pot_dy   += source_cluster_q_dy[j]   * r3inv * dy;
            temp_pot_dz   += source_cluster_q_dz[j]   * r3inv * dz;
            temp_pot_dxy  += source_cluster_q_dxy[j]  * r5inv * dx * dy;
            temp_pot_dyz  += source_cluster_q_dyz[j]  * r5inv * dy * dz;
            temp_pot_dxz  += source_cluster_q_dxz[j]  * r5inv * dx * dz;
            temp_pot_dxyz += source_cluster_q_dxyz[j] * r7inv * dx * dy * dz * 15.0;

        } // end loop over interpolation points
        
#ifdef OPENACC_ENABLED
        #pragma acc atomic
        target_cluster_q_[i]     += temp_pot_;
        #pragma acc atomic
        target_cluster_q_dx[i]   += temp_pot_dx;
        #pragma acc atomic
        target_cluster_q_dy[i]   += temp_pot_dy;
        #pragma acc atomic
        target_cluster_q_dz[i]   += temp_pot_dz;
        #pragma acc atomic
        target_cluster_q_dxy[i]  += temp_pot_dxy;
        #pragma acc atomic
        target_cluster_q_dyz[i]  += temp_pot_dyz;
        #pragma acc atomic
        target_cluster_q_dxz[i]  += temp_pot_dxz;
        #pragma acc atomic
        target_cluster_q_dxyz[i] += temp_pot_dxyz;
#else
        target_cluster_q_[i]     += temp_pot_;
        target_cluster_q_dx[i]   += temp_pot_dx;
        target_cluster_q_dy[i]   += temp_pot_dy;
        target_cluster_q_dz[i]   += temp_pot_dz;
        target_cluster_q_dxy[i]  += temp_pot_dxy;
        target_cluster_q_dyz[i]  += temp_pot_dyz;
        target_cluster_q_dxz[i]  += temp_pot_dxz;
        target_cluster_q_dxyz[i] += temp_pot_dxyz;
#endif


    }
#ifdef OPENACC_ENABLED
    } // end kernel
#endif
    return;
}
