/* Interaction Kernels */
#ifndef H_DCF_H
#define H_DCF_H
 

void dcfDirect(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster,
        double *target_x, double *target_y, double *target_z,
        double *source_x, double *source_y, double *source_z, double *source_q, double *source_w,
        double kernel_parameter, double *potential, int gpu_async_stream_id);

void dcfApproximationLagrange(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster,
        double *target_x, double *target_y, double *target_z,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
        double kernel_parameter, double *potential, int gpu_async_stream_id);

void dcfApproximationHermite(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster, int total_number_interpolation_points,
        double *target_x, double *target_y, double *target_z,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
        double kernel_parameter, double *potential, int gpu_async_stream_id);

#endif /* H_DCF_H */
