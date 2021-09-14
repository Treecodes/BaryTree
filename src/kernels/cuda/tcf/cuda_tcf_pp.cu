#include <math.h>
#include <float.h>
#include <stdio.h>

#include "cuda_tcf_pp.h"


__global__ 
static void compute()
{
    printf("Hello TCF_PP thread %d, block %d\n", threadIdx.x, blockIdx.x);
}


__host__
void K_CUDA_TCF_PP(
    int target_x_low_ind,  int target_x_high_ind,
    int target_y_low_ind,  int target_y_high_ind,
    int target_z_low_ind,  int target_z_high_ind,
    double target_xmin,    double target_ymin,    double target_zmin,
    double target_xdd,     double target_ydd,     double target_zdd,
    int target_x_dim_glob, int target_y_dim_glob, int target_z_dim_glob,
    int cluster_num_sources, int cluster_idx_start,
    double *source_x, double *source_y, double *source_z, double *source_q,
    struct RunParams *run_params, double *potential, int gpu_async_stream_id)
{
    double kap = run_params->kernel_params[0];
    double eta = run_params->kernel_params[1];
    double kap_eta_2 = kap * eta / 2.0;

    compute<<<1,32>>>();
    cudaDeviceSynchronize();

    return;
}
