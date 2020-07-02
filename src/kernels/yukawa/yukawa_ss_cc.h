/* Interaction Kernels */
#ifndef H_K_YUKAWA_SS_CC_H
#define H_K_YUKAWA_SS_CC_H
 
#include "../../run_params/struct_run_params.h"


void K_Yukawa_SS_CC_Lagrange(int number_of_sources_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_sources, int starting_index_of_cluster,
        double *source_cluster_x, double *source_cluster_y, double *source_cluster_z, double *source_cluster_q, double *source_cluster_w,
        double *target_cluster_x, double *target_cluster_y, double *target_cluster_z, double *target_cluster_charge, double *target_cluster_w,
        struct RunParams *run_params, int gpu_async_stream_id);


#endif /* H_K_YUKAWA_SS_CC_H */
