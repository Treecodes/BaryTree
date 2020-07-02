/* Interaction Kernels */
#ifndef H_K_REGULARIZED_COULOMB_SS_PP_H
#define H_K_REGULARIZED_COULOMB_SS_PP_H
 
#include "../../run_params/struct_run_params.h"


void K_RegularizedCoulomb_SS_PP(int number_of_targets_in_batch, int number_of_source_points_in_cluster,
        int starting_index_of_target, int starting_index_of_source,
        double *target_x, double *target_y, double *target_z, double *target_charge,
        double *source_x, double *source_y, double *source_z, double *source_charge, double * source_weight,
        struct RunParams *run_params, double *potential, int gpu_async_stream_id);


#endif /* H_K_REGULARIZED_COULOMB_SS_PP_H */
