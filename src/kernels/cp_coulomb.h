/* Interaction Kernels */
#ifndef H_COULOMB_H
#define H_COULOMB_H
 
#include "../struct_kernel.h"


void CP_coulombApproximationLagrange(int number_of_sources_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_sources, int starting_index_of_cluster,
        double *source_x, double *source_y, double *source_z, double *source_q,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
        struct kernel *kernel, int gpu_async_stream_id);

void CP_coulombApproximationHermite(int number_of_sources_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_sources, int starting_index_of_cluster, int total_number_interpolation_points,
        double *source_x, double *source_y, double *source_z, double *source_q,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
        struct kernel *kernel, int gpu_async_stream_id);

#endif /* H_COULOMB_H */
