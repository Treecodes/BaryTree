#include <math.h>
#include <float.h>
#include <stdio.h>

void coulombSingularitySubtractionDirect( int number_of_targets_in_batch, int number_of_source_points_in_cluster, int starting_index_of_target, int starting_index_of_source,
                    double *target_x, double *target_y, double *target_z, double *target_charge,
                    double *source_x, double *source_y, double *source_z, double *source_charge, double * source_weight,
                    double kernel_parameter, double *potential, int gpu_async_stream_id){

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id)
    #pragma acc loop independent
    {
#endif
    for (int i = 0; i < number_of_targets_in_batch; i++) {
        int ii=starting_index_of_target + i;
        double temporary_potential = 0.0;

        for (int j = 0; j < number_of_source_points_in_cluster; j++) {
            int jj=starting_index_of_source + j;
            double dx = target_x[ ii] - source_x[ jj];
            double dy = target_y[ ii] - source_y[ jj];
            double dz = target_z[ ii] - source_z[ jj];
            double r  = sqrt( dx*dx + dy*dy + dz*dz);

            if (r > DBL_MIN){
                temporary_potential += ( source_charge[jj] - target_charge[ii] * exp(-r*r/kernel_parameter/kernel_parameter) ) * source_weight[jj] / r;
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


void coulombSingularitySubtractionApproximationLagrange( int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster, int starting_index_of_target, int starting_index_of_cluster,
                    double *target_x, double *target_y, double *target_z, double *target_charge,
                    double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge, double *cluster_weight,
                    double kernel_parameter, double *potential, int gpu_async_stream_id){


#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id)
    #pragma acc loop independent
    {
#endif
    for (int i = 0; i < number_of_targets_in_batch; i++) {
        int ii=starting_index_of_target + i;
        double temporary_potential = 0.0;

        for (int j = 0; j < number_of_interpolation_points_in_cluster; j++) {
            int jj=starting_index_of_cluster + j;
            double dx = target_x[ ii] - cluster_x[ jj];
            double dy = target_y[ ii] - cluster_y[ jj];
            double dz = target_z[ ii] - cluster_z[ jj];
            double r  = sqrt( dx*dx + dy*dy + dz*dz);

            if (r > DBL_MIN){
                temporary_potential += ( cluster_charge[jj] - target_charge[ii]*cluster_weight[jj] * exp(-r*r/kernel_parameter/kernel_parameter) ) / r;
//                        return (sourceQ - targetQ * sourceW) * G;

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
