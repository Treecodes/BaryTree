#include <math.h>
#include <float.h>
#include <stdio.h>

void yukawaDirect( int number_of_targets_in_batch, int number_of_source_points_in_cluster, int starting_index_of_target, int starting_index_of_source,
                    double *target_x, double *target_y, double *target_z,
                    double *source_x, double *source_y, double *source_z, double *source_charge, double * source_weight,
                    double kernel_parameter, double *potential, int gpu_async_stream_id){

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id)
    #pragma acc loop independent
    {
#endif
    for (int i = 0; i < number_of_targets_in_batch; i++) {

        double temporary_potential = 0.0;

        for (int j = 0; j < number_of_source_points_in_cluster; j++) {

            double dx = target_x[ starting_index_of_target + i] - source_x[ starting_index_of_source + j];
            double dy = target_y[ starting_index_of_target + i] - source_y[ starting_index_of_source + j];
            double dz = target_z[ starting_index_of_target + i] - source_z[ starting_index_of_source + j];
            double r  = sqrt( dx*dx + dy*dy + dz*dz);

            if (r > DBL_MIN){
                temporary_potential += source_charge[starting_index_of_source + j] * source_weight[starting_index_of_source + j] * exp(-kernel_parameter*r) / r;
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


void yukawaApproximationLagrange( int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster, int starting_index_of_target, int starting_index_of_cluster,
                    double *target_x, double *target_y, double *target_z,
                    double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
                    double kernel_parameter, double *potential, int gpu_async_stream_id){


#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id)
    #pragma acc loop independent
    {
#endif
    for (int i = 0; i < number_of_targets_in_batch; i++) {

        double temporary_potential = 0.0;

        for (int j = 0; j < number_of_interpolation_points_in_cluster; j++) {

            double dx = target_x[ starting_index_of_target + i] - cluster_x[ starting_index_of_cluster + j];
            double dy = target_y[ starting_index_of_target + i] - cluster_y[ starting_index_of_cluster + j];
            double dz = target_z[ starting_index_of_target + i] - cluster_z[ starting_index_of_cluster + j];
            double r  = sqrt( dx*dx + dy*dy + dz*dz);

            if (r > DBL_MIN){
                temporary_potential += cluster_charge[starting_index_of_cluster + j] * exp(-kernel_parameter*r) /r;
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
