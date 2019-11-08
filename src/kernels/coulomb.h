/* Interaction Kernels */
#ifndef H_COULOMB_H
#define H_COULOMB_H
 

void coulombDirect(                 int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster, int starting_index_of_target, int starting_index_of_cluster,
                                    double *target_x, double *target_y, double *target_z,
                                    double *source_x, double *source_y, double *source_z, double *source_charge, double * source_weight,
                                    double *potential, int gpu_async_stream_id);

void coulombApproximationLagrange(  int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster, int starting_index_of_target, int starting_index_of_cluster,
                                    double *target_x, double *target_y, double *target_z,
                                    double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
                                    double *potential, int gpu_async_stream_id);

void coulombApproximationHermite(   int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster, int starting_index_of_target,
                                    int starting_index_of_cluster, int total_number_interpolation_points,
                                    double *target_x, double *target_y, double *target_z,
                                    double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
                                    double *potential, int gpu_async_stream_id);

#endif /* H_COULOMB_H */
