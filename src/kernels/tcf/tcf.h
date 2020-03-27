/* Interaction Kernels */
#ifndef H_TCF_H
#define H_TCF_H
 

void tcfDirect(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster,
        double *target_x, double *target_y, double *target_z,
        double *source_x, double *source_y, double *source_z, double *source_q, double *source_w,
        double kernel_parameter1, double kernel_parameter2, double *potential, int gpu_async_stream_id);

void tcfApproximationLagrange(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster,
        double *target_x, double *target_y, double *target_z,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
        double kernel_parameter1, double kernel_parameter2, double *potential, int gpu_async_stream_id);

#endif /* H_TCF_H */
