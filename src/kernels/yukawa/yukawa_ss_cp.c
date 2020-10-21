#include <math.h>
#include <float.h>
#include <stdio.h>

#include "../../run_params/struct_run_params.h"
#include "yukawa_ss_cp.h"


void K_Yukawa_SS_CP_Lagrange(int number_of_sources_in_batch, int number_of_interpolation_points_in_cluster,
         int starting_index_of_sources, int starting_index_of_cluster,
         double *source_x, double *source_y, double *source_z, double *source_q, double *source_w,
         double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q, double *cluster_w,
         struct RunParams *run_params, int gpu_async_stream_id)
{

    double kernel_parameter = run_params->kernel_params[0];

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(source_x, source_y, source_z, source_q, source_w,\
                        cluster_x, cluster_y, cluster_z, cluster_q, cluster_w)
    {
#endif
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif	
    for (int i = 0; i < number_of_interpolation_points_in_cluster; i++) {

        double temporary_potential = 0.0;
        double temporary_weight = 0.0;

        double cx = cluster_x[starting_index_of_cluster + i];
        double cy = cluster_y[starting_index_of_cluster + i];
        double cz = cluster_z[starting_index_of_cluster + i];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential) reduction(+:temporary_weight)
#endif
        for (int j = 0; j < number_of_sources_in_batch; j++) {
#ifdef OPENACC_ENABLED
            #pragma acc cache(source_x[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_y[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_z[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_q[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch], \
                              source_w[starting_index_of_sources : starting_index_of_sources+number_of_sources_in_batch])
#endif

            int jj = starting_index_of_sources + j;
            double dx = cx - source_x[jj];
            double dy = cy - source_y[jj];
            double dz = cz - source_z[jj];
            double r = sqrt(dx*dx + dy*dy + dz*dz);

            temporary_potential += source_q[jj]  * source_w[jj] * exp(-kernel_parameter*r) /r;
            temporary_weight    +=                 source_w[jj] * exp(-kernel_parameter*r) /r;

        } // end loop over interpolation points
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        cluster_q[starting_index_of_cluster + i] += temporary_potential;
#ifdef OPENACC_ENABLED
        #pragma acc atomic
#endif
        cluster_w[starting_index_of_cluster + i] += temporary_weight;
    }
#ifdef OPENACC_ENABLED
    } // end kernel
#endif
    return;
}


