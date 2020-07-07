/* Interaction Kernels */
#ifndef H_K_TCF_DIRECT_H
#define H_K_TCF_DIRECT_H
 
#include "../../run_params/struct_run_params.h"


void K_TCF_Direct(int target_x_low_ind,  int target_x_high_ind,
                  int target_y_low_ind,  int target_y_high_ind,
                  int target_z_low_ind,  int target_z_high_ind,
                  double target_xmin,    double target_ymin,    double target_zmin,
                  
                  double target_xdd,     double target_ydd,     double target_zdd,
                  int target_x_dim_glob, int target_y_dim_glob, int target_z_dim_glob,

                  int number_of_source_points_in_cluster, int starting_index_of_source,
                  double *source_x, double *source_y, double *source_z, double *source_q,

                  struct RunParams *run_params, double *potential, int gpu_async_stream_id);

//void K_TCF_Direct(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
//        int starting_index_of_target, int starting_index_of_cluster,
//        double *target_x, double *target_y, double *target_z,
//        double *source_x, double *source_y, double *source_z, double *source_charge, double *source_weight,
//        struct RunParams *run_params, double *potential, int gpu_async_stream_id);


#endif /* H_K_TCF_DIRECT_H */
