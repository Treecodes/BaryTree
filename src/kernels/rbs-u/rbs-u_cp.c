#include <math.h>
#include <float.h>
#include <stdio.h>

#include "../../run_params/struct_run_params.h"
#include "rbs-u_cp.h"


void K_RBSu_CP_Lagrange(int number_of_sources_in_batch, int number_of_interpolation_points_in_cluster,
         int starting_index_of_source, int starting_index_of_cluster,
         double *source_x, double *source_y, double *source_z, double *source_q,
         double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
         struct RunParams *run_params, int gpu_async_stream_id)
{

    double delta = run_params->kernel_params[0];
    double delta2 = delta * delta;

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

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential)
#endif
        for (int j = 0; j < number_of_sources_in_batch; j++) {
#ifdef OPENACC_ENABLED
            #pragma acc cache(source_x[starting_index_of_source : starting_index_of_source+number_of_sources_in_batch], \
                              source_y[starting_index_of_source : starting_index_of_source+number_of_sources_in_batch], \
                              source_z[starting_index_of_source : starting_index_of_source+number_of_sources_in_batch], \
                              source_q[starting_index_of_source : starting_index_of_source+number_of_sources_in_batch])
#endif

            int jj = starting_index_of_source + j;
            double dx = cx - source_x[jj];
            double dy = cy - source_y[jj];
            double r = dx*dx + dy*dy + delta2;

            temporary_potential += 1. / (2. * M_PI) * source_q[jj] * dy / r;

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
