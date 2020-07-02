/* Interaction Kernels */
#ifndef H_K_REGULARIZED_YUKAWA_SS_PC_H
#define H_K_REGULARIZED_YUKAWA_SS_PC_H
 
#include "../../run_params/struct_run_params.h"


void K_RegularizedYukawa_SS_PC_Lagrange(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster,
        double *target_x, double *target_y, double *target_z, double *target_charge, double *cluster_weight,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
        struct RunParams *run_params, double *potential, int gpu_async_stream_id);



#endif /* H_K_REGULARIZED_YUKAWA_SS_PC_H */
