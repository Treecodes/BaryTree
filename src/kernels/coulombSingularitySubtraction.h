/* Interaction Kernels */
#ifndef H_COULOMBSINGULARITYSUBTRACTION_H
#define H_COULOMBSINGULARITYSUBTRACTION_H
 

void coulombSingularitySubtractionDirect(int number_of_targets_in_batch, int number_of_source_points_in_cluster, int starting_index_of_target, int starting_index_of_source,
                                        double *target_x, double *target_y, double *target_z, double *target_charge,
                                        double *source_x, double *source_y, double *source_z, double *source_charge, double * source_weight,
                                        double kernel_parameter, double *potential, int gpu_async_stream_id);

void coulombSingularitySubtractionApproximationLagrange( int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster, int starting_index_of_target, int starting_index_of_cluster,
                                                        double *target_x, double *target_y, double *target_z, double *target_charge, double *cluster_weight,
                                                        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
                                                        double kernel_parameter, double *potential, int gpu_async_stream_id);

#endif /* H_COULOMBSINGULARITYSUBTRACTION_H */
