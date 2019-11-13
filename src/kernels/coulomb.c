#include <math.h>
#include <float.h>
#include <stdio.h>

void coulombDirect(int number_of_targets_in_batch, int number_of_source_points_in_cluster,
                   int starting_index_of_target, int starting_index_of_source,
                   double *target_x, double *target_y, double *target_z,
                   double *source_x, double *source_y, double *source_z, double *source_charge, double *source_weight,
                   double *potential, int gpu_async_stream_id)
{

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x,target_y,target_z,source_x,source_y,source_z,source_charge,source_weight,potential)
    {
    #pragma acc loop independent
#endif
    for (int i = 0; i < number_of_targets_in_batch; i++) {
        int ii=starting_index_of_target + i;
        double temporary_potential = 0.0;

        double tx = target_x[ii];
        double ty = target_y[ii];
        double tz = target_z[ii];

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential)
#endif
        for (int j = 0; j < number_of_source_points_in_cluster; j++) {
            int jj=starting_index_of_source + j;
            double dx = tx - source_x[jj];
            double dy = ty - source_y[jj];
            double dz = tz - source_z[jj];
            double r  = sqrt(dx*dx + dy*dy + dz*dz);

            if (r > DBL_MIN) {
                temporary_potential += source_charge[jj] * source_weight[jj] / r;
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


void coulombApproximationLagrange(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
                                  int starting_index_of_target, int starting_index_of_cluster,
                                  double *target_x, double *target_y, double *target_z,
                                  double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
                                  double *potential, int gpu_async_stream_id)
{

#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x,target_y,target_z,cluster_x,cluster_y,cluster_z,cluster_charge,potential)
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

            double dx = tx - cluster_x[ starting_index_of_cluster + j];
            double dy = ty - cluster_y[ starting_index_of_cluster + j];
            double dz = tz - cluster_z[ starting_index_of_cluster + j];
            double r  = sqrt(dx*dx + dy*dy + dz*dz);

            if (r > DBL_MIN){
                temporary_potential += cluster_charge[starting_index_of_cluster + j]/r;
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

void coulombApproximationHermite( int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster, int starting_index_of_target,
                    int starting_index_of_cluster, int total_number_interpolation_points,
                    double *target_x, double *target_y, double *target_z,
                    double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
                    double *potential, int gpu_async_stream_id){


    // total_number_interpolation_points is the stride, separating clustersQ, clustersQx, clustersQy, etc.
    double *cluster_delta_x =   &cluster_charge[1*total_number_interpolation_points];
    double *cluster_delta_y =   &cluster_charge[2*total_number_interpolation_points];
    double *cluster_delta_z =   &cluster_charge[3*total_number_interpolation_points];
    double *cluster_delta_xy =  &cluster_charge[4*total_number_interpolation_points];
    double *cluster_delta_yz =  &cluster_charge[5*total_number_interpolation_points];
    double *cluster_delta_xz =  &cluster_charge[6*total_number_interpolation_points];
    double *cluster_delta_xyz = &cluster_charge[7*total_number_interpolation_points];


#ifdef OPENACC_ENABLED
    #pragma acc kernels async(gpu_async_stream_id) present(target_x,target_y,target_z,cluster_x,cluster_y,cluster_z,cluster_charge,potential, \
            cluster_delta_x,cluster_delta_y,cluster_delta_z,cluster_delta_xy,cluster_delta_yz,cluster_delta_xz,cluster_delta_xyz)
    {
#endif
#ifdef OPENACC_ENABLED
    #pragma acc loop independent
#endif
    for (int i = 0; i < number_of_targets_in_batch; i++) {

        int ii=starting_index_of_target + i;
        double temporary_potential = 0.0;

#ifdef OPENACC_ENABLED
        #pragma acc loop independent reduction(+:temporary_potential)
#endif
        for (int j = 0; j < number_of_interpolation_points_in_cluster; j++) {
            int jj=starting_index_of_cluster + j;
            double dx = target_x[ii] - cluster_x[jj];
            double dy = target_y[ii] - cluster_y[jj];
            double dz = target_z[ii] - cluster_z[jj];
            double r  = sqrt( dx*dx + dy*dy + dz*dz);

            double rinv  = 1 / r;
            double r3inv = rinv*rinv*rinv;
            double r5inv = r3inv*rinv*rinv;
            double r7inv = r5inv*rinv*rinv;

            if (r > DBL_MIN){

                temporary_potential +=  rinv  * ( cluster_charge[jj])
                                +      r3inv * ( cluster_delta_x[jj]*dx +  cluster_delta_y[jj]*dy +  cluster_delta_z[jj]*dz )
                                + 3 *  r5inv * ( cluster_delta_xy[jj]*dx*dy +  cluster_delta_yz[jj]*dy*dz +  cluster_delta_xz[jj]*dx*dz )
                                + 15 * r7inv *   cluster_delta_xyz[jj]*dx*dy*dz;

//                temporary_potential += cluster_charge[starting_index_of_cluster + j]/r;
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