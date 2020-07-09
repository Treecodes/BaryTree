/* Interaction Kernels */
#ifndef H_K_TCF_CP_H
#define H_K_TCF_CP_H
 
#include "../../run_params/struct_run_params.h"


void K_TCF_CP_Lagrange(int number_of_sources_in_batch, int starting_index_of_sources,
                       int interp_pts_per_cluster, int cluster_charge_start,
                       int interp_order_lim, int cluster_pts_start,
         double *source_x, double *source_y, double *source_z, double *source_q,
         double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
         struct RunParams *run_params, int gpu_async_stream_id);

void K_TCF_CP_Hermite(int number_of_sources_in_batch, int number_of_interpolation_points_in_cluster,
        int starting_index_of_sources, int starting_index_of_cluster,
        double *source_x, double *source_y, double *source_z, double *source_q,
        double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_charge,
        struct RunParams *run_params, int gpu_async_stream_id);


#endif /* H_K_TCF_CP_H */
