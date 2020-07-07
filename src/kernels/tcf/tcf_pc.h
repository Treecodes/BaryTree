/* Interaction Kernels */
#ifndef H_K_TCF_PC_H
#define H_K_TCF_PC_H
 
#include "../../run_params/struct_run_params.h"


//void K_TCF_PC_Lagrange(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
//        int starting_index_of_target, int starting_index_of_cluster,
//        double *target_x, double *target_y, double *target_z,
//        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
//        struct RunParams *run_params, double *potential, int gpu_async_stream_id);

void K_TCF_PC_Lagrange(int target_x_low_ind,  int target_x_high_ind,
                       int target_y_low_ind,  int target_y_high_ind,
                       int target_z_low_ind,  int target_z_high_ind,
                       double target_xmin,    double target_ymin,    double target_zmin,
                       
                       double target_xdd,     double target_ydd,     double target_zdd,
                       int target_x_dim_glob, int target_y_dim_glob, int target_z_dim_glob,

                       int number_of_interpolation_points_in_cluster, int starting_index_of_cluster,
                       double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,

                       struct RunParams *run_params, double *potential, int gpu_async_stream_id);

//void K_TCF_PC_Hermite(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
//        int starting_index_of_target, int starting_index_of_cluster, int total_number_interpolation_points,
//        double *target_x, double *target_y, double *target_z,
//        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
//        struct RunParams *run_params, double *potential, int gpu_async_stream_id);

void K_TCF_PC_Hermite(int target_x_low_ind,  int target_x_high_ind,
                      int target_y_low_ind,  int target_y_high_ind,
                      int target_z_low_ind,  int target_z_high_ind,
                      double target_xmin,    double target_ymin,    double target_zmin,
                      
                      double target_xdd,     double target_ydd,     double target_zdd,
                      int target_x_dim_glob, int target_y_dim_glob, int target_z_dim_glob,

                      int number_of_interpolation_points_in_cluster, int starting_index_of_cluster,
                      double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,

                      struct RunParams *run_params, double *potential, int gpu_async_stream_id);


#endif /* H_K_TCF_PC_H */
