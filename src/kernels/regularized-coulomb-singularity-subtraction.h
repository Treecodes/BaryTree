/* Interaction Kernels */
#ifndef H_REGULARIZEDCOULOMBSINGULARITYSUBTRACTION_H
#define H_REGULARIZEDCOULOMBSINGULARITYSUBTRACTION_H
 
#include "../struct_kernel.h"


void regularizedCoulombSingularitySubtractionDirect(int number_of_targets_in_batch, int number_of_source_points_in_cluster,
        int starting_index_of_target, int starting_index_of_source,
        double *target_x, double *target_y, double *target_z, double *target_charge,
        double *source_x, double *source_y, double *source_z, double *source_charge, double * source_weight,
        struct kernel *kernel, double *potential, int gpu_async_stream_id);

void regularizedCoulombSingularitySubtractionApproximationLagrange(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster,
        double *target_x, double *target_y, double *target_z, double *target_charge, double *cluster_weight,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
        struct kernel *kernel, double *potential, int gpu_async_stream_id);

void regularizedCoulombSingularitySubtractionApproximationHermite(int number_of_targets_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_target, int starting_index_of_cluster, int total_number_interpolation_points,
        double *target_x, double *target_y, double *target_z, double *target_charge,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge, double *cluster_weight,
        struct kernel *kernel, double *potential, int gpu_async_stream_id);

void regularizedCoulombSingularitySubtractionCorrection(double *potential, double *target_q,
        int numTargets, struct kernel *kernel);


#endif /* H_REGULARIZEDCOULOMBSINGULARITYSUBTRACTION_H */
