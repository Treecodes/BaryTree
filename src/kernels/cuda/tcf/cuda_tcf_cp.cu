#include <math.h>
#include <float.h>
#include <stdio.h>

#include "cuda_tcf_cp.h"


__global__ 
static void compute(double kap, int batch_num_sources, int batch_idx_start,
    int cluster_q_start, int cluster_pts_start, int interp_order_lim,
    double *source_x, double *source_y, double *source_z, double *source_q,
    double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q)
{
    printf("Hello TCF_CP thread %d, block %d\n", threadIdx.x, blockIdx.x);
}


__host__
void K_CUDA_TCF_CP_Lagrange(
    int batch_num_sources, int batch_idx_start, 
    int cluster_q_start, int cluster_pts_start, int interp_order_lim,
    double *source_x, double *source_y, double *source_z, double *source_q,
    double *cluster_x, double *cluster_y, double *cluster_z, double *cluster_q,
    struct RunParams *run_params, int gpu_async_stream_id)
{
    double kap = run_params->kernel_params[0];
    //double eta = run_params->kernel_params[1];
    //double kap_eta_2 = kap * eta / 2.0;

    compute<<<1,32>>>(kap, batch_num_sources, batch_idx_start,
                      cluster_q_start, cluster_pts_start, interp_order_lim,
                      source_x,  source_y,  source_z,  source_q,
                      cluster_x, cluster_y, cluster_z, cluster_q);
    cudaDeviceSynchronize();

    return;
}
