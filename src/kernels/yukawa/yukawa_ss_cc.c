#include <math.h>
#include <float.h>
#include <stdio.h>

#include "../../run_params/struct_run_params.h"
#include "yukawa_ss_cc.h"


void K_Yukawa_SS_CC_Lagrange(int number_of_sources_in_batch, int number_of_interpolation_points_in_cluster,
         int starting_index_of_sources, int starting_index_of_cluster,
         double *source_cluster_x, double *source_cluster_y, double *source_cluster_z, double *source_cluster_q, double *source_cluster_w,
         double *target_cluster_x, double *target_cluster_y, double *target_cluster_z, double *target_cluster_q, double *target_cluster_w,
         struct RunParams *run_params, int gpu_async_stream_id)
{

    double kernel_parameter = run_params->kernel_params[0];

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(source_cluster_x, source_cluster_y, source_cluster_z, source_cluster_q, source_cluster_w, \
                        target_cluster_x, target_cluster_y, target_cluster_z, target_cluster_q, target_cluster_w)
    {
#endif
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif	
    for (int i = 0; i < number_of_interpolation_points_in_cluster; i++) {

        double temporary_potential = 0.0;
        double temporary_weight = 0.0;

        double cx = target_cluster_x[starting_index_of_cluster + i];
        double cy = target_cluster_y[starting_index_of_cluster + i];
        double cz = target_cluster_z[starting_index_of_cluster + i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential,temporary_weight)
#endif
        for (int j = 0; j < number_of_sources_in_batch; j++) {
#ifdef OPENACC_ENABLED
            #pragma acc cache(source_cluster_x[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_cluster_y[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_cluster_z[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_cluster_q[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_cluster_w[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch])
#endif

            int jj = starting_index_of_sources + j;
            double dx = cx - source_cluster_x[jj];
            double dy = cy - source_cluster_y[jj];
            double dz = cz - source_cluster_z[jj];
            double r = sqrt(dx*dx + dy*dy + dz*dz);

            temporary_potential += source_cluster_q[jj] * exp(-kernel_parameter*r) /r; // source_cluster_q already has source_q * source_w
            temporary_weight    += source_cluster_w[jj] * exp(-kernel_parameter*r) /r;

        } // end loop over interpolation points
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        target_cluster_q[starting_index_of_cluster + i] += temporary_potential;
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        target_cluster_w[starting_index_of_cluster + i] += temporary_weight;
    }
#ifdef OPENACC_ENABLED
    } // end kernel
#endif
    return;
}


