#include <math.h>
#include <float.h>
#include <stdio.h>

void yukawaSingularitySubtractionDirect( int number_of_targets_in_batch, int number_of_source_points_in_cluster, int starting_index_of_target, int starting_index_of_source,
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
                temporary_potential += ( source_charge[jj] - target_charge[ii]) * source_weight[jj] * exp(-kernel_parameter*r) / r;
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


void yukawaSingularitySubtractionApproximationLagrange( int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster, int starting_index_of_target, int starting_index_of_cluster,
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
                temporary_potential += ( cluster_charge[jj] - target_charge[ii]*cluster_weight[jj] ) * exp(-kernel_parameter*r) /r;
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


void yukawaSingularitySubtractionApproximationHermite( int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster, int starting_index_of_target,
                    int starting_index_of_cluster, int total_number_interpolation_points,
                    double *target_x, double *target_y, double *target_z, double *target_charge,
                    double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge, double *cluster_weight,
                    double kernel_parameter, double *potential, int gpu_async_stream_id){


    // total_number_interpolation_points is the stride, separating clustersQ, clustersQx, clustersQy, etc.
    double *cluster_charge_delta_x =   &cluster_charge[1*total_number_interpolation_points];
    double *cluster_charge_delta_y =   &cluster_charge[2*total_number_interpolation_points];
    double *cluster_charge_delta_z =   &cluster_charge[3*total_number_interpolation_points];
    double *cluster_charge_delta_xy =  &cluster_charge[4*total_number_interpolation_points];
    double *cluster_charge_delta_yz =  &cluster_charge[5*total_number_interpolation_points];
    double *cluster_charge_delta_xz =  &cluster_charge[6*total_number_interpolation_points];
    double *cluster_charge_delta_xyz = &cluster_charge[7*total_number_interpolation_points];

    double *cluster_weight_delta_x =   &cluster_weight[1*total_number_interpolation_points];
    double *cluster_weight_delta_y =   &cluster_weight[2*total_number_interpolation_points];
    double *cluster_weight_delta_z =   &cluster_weight[3*total_number_interpolation_points];
    double *cluster_weight_delta_xy =  &cluster_weight[4*total_number_interpolation_points];
    double *cluster_weight_delta_yz =  &cluster_weight[5*total_number_interpolation_points];
    double *cluster_weight_delta_xz =  &cluster_weight[6*total_number_interpolation_points];
    double *cluster_weight_delta_xyz = &cluster_weight[7*total_number_interpolation_points];


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
            double dx = target_x[ii] - cluster_x[jj];
            double dy = target_y[ii] - cluster_y[jj];
            double dz = target_z[ii] - cluster_z[jj];
            double r  = sqrt( dx*dx + dy*dy + dz*dz);

            double r2 = r*r;
            double r3 = r2*r;
            double r4 = r2*r2;
            double rinv = 1 / r;
            double r3inv = rinv*rinv*rinv;
            double r5inv = r3inv*rinv*rinv;
            double r7inv = r5inv*rinv*rinv;
            double kernel_parameter2=kernel_parameter*kernel_parameter;
            double kernel_parameter3=kernel_parameter*kernel_parameter*kernel_parameter;

            double charge_diff = cluster_charge[jj] - cluster_weight[ii];
            double delta_x_diff = cluster_charge_delta_x[jj]-cluster_weight_delta_x[ii];
            double delta_y_diff = cluster_charge_delta_y[jj]-cluster_weight_delta_y[ii];
            double delta_z_diff = cluster_charge_delta_z[jj]-cluster_weight_delta_z[ii];
            double delta_xy_diff = cluster_charge_delta_xy[jj]-cluster_weight_delta_xy[ii];
            double delta_yz_diff = cluster_charge_delta_yz[jj]-cluster_weight_delta_yz[ii];
            double delta_xz_diff = cluster_charge_delta_xz[jj]-cluster_weight_delta_xz[ii];
            double delta_xyz_diff = cluster_charge_delta_xyz[jj]-cluster_weight_delta_xyz[ii];

            if (r > DBL_MIN){

                temporary_potential +=       exp(-kernel_parameter*r)*(
                                   rinv  * ( charge_diff)
                            +      r3inv * (1 + kernel_parameter*r) * ( delta_x_diff*dx +  delta_y_diff*dy +  delta_z_diff*dz )
                            +      r5inv * (3 + 3*kernel_parameter*r + kernel_parameter2*r2 ) * ( delta_xy_diff*dx*dy +  delta_yz_diff*dy*dz +  delta_xz_diff*dx*dz )
                            +      r7inv * (15 + 15*kernel_parameter*r + 6*kernel_parameter2*r2 + kernel_parameter3*r3) * delta_xyz_diff*dx*dy*dz);




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
