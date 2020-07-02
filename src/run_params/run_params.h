#ifndef H_RUN_PARAMS_FUNCTIONS_H
#define H_RUN_PARAMS_FUNCTIONS_H

#include "../utilities/enums.h"
#include "struct_run_params.h"


void RunParams_Setup(struct RunParams **run_params_addr,
                     KERNEL kernel, int num_kernel_params, double *kernel_params,
                     APPROXIMATION approximation,
                     SINGULARITY singularity,
                     COMPUTE_TYPE compute_type,
                     double theta, int interp_degree,
                     int max_per_source_leaf, int max_per_target_leaf, double size_check_factor,
                     double beta, int verbosity);

void RunParams_Validate(struct RunParams *run_params);

void RunParams_Free(struct RunParams **run_params_addr);

void RunParams_Print(struct RunParams *run_params);


#endif
