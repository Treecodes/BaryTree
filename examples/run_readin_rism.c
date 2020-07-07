#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <mpi.h>

#include "../src/utilities/tools.h"
#include "../src/utilities/timers.h"

#include "../src/particles/struct_particles.h"
#include "../src/run_params/struct_run_params.h"
#include "../src/run_params/run_params.h"

#include "../src/drivers/treedriver.h"
#include "../src/drivers/directdriver.h"

#include "support_fns.h"


int main(int argc, char **argv)
{
    /* MPI initialization */
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    if (rank == 0) printf("[random cube example] Beginning random cube example with %d ranks.\n", numProcs);

    /* run parameters */
    int N, run_direct, slice;
    double xyz_limits[6];
    int grid_dim[3];
    double grid_dd[3];
    FILE *fp_pqr;
    char file_pqr[256];
    
    struct RunParams *run_params = NULL;
    
    FILE *fp = fopen(argv[1], "r");
    Params_Parse_Readin(fp, &run_params, &N, file_pqr, &run_direct, &slice, xyz_limits, grid_dim);
    
    //Slicing currently doesn't work! So I set it to 1.
    slice = 1;

    grid_dd[0] = (xyz_limits[1] - xyz_limits[0]) / (grid_dim[0] - 1);
    grid_dd[1] = (xyz_limits[3] - xyz_limits[2]) / (grid_dim[1] - 1);
    grid_dd[2] = (xyz_limits[5] - xyz_limits[4]) / (grid_dim[2] - 1);


    /* data structures for BaryTree calculation and comparison */
    struct Particles *sources = NULL;
    struct Particles *targets = NULL;
    //struct Particles *targets_sample = NULL;
    double *potential = NULL, *potential_direct = NULL;
    
    /* variables for collecting accuracy info */
    double potential_engy = 0, potential_engy_glob = 0;
    double potential_engy_direct = 0, potential_engy_direct_glob = 0;
    double glob_inf_err = 0, glob_n2_err = 0, glob_relinf_err = 0, glob_reln2_err = 0;

    /* variables for date-time calculation */
    double time_run[4], time_tree[13], time_direct[4];
    double time_run_glob[3][4], time_tree_glob[3][13], time_direct_glob[3][4];


    /* Beginning total runtime timer */
    START_TIMER(&time_run[3]);
    
    
    //~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Setup
    //~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    START_TIMER(&time_run[0]);


    /* Setting up sources with MPI-allocated source arrays for RMA use */
    
    sources = malloc(sizeof(struct Particles));
    sources->num = N;

    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->x));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->y));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->z));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->q));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->w));

    char c[256], c1[120], c2[120], c3[120], c4[10], c5[10];
    double a1, a2, a3, b1, b2;
    int atom_ctr = 0;

    FILE *points_fp = fopen(file_pqr, "r");

//    while (fgets(c, sizeof(c), points_fp)) {
//        sscanf(c, "%s %s %s %s %s %lf %lf %lf %lf %lf",
//               c1, c2, c3, c4, c5, &a1, &a2, &a3, &b1, &b2);
//
//        if (strncmp(c1, "ATOM", 4) == 0) {
//            sources->x[atom_ctr] = a1;
//            sources->y[atom_ctr] = a2;
//            sources->z[atom_ctr] = a3;
//            sources->q[atom_ctr] = b2;
//            sources->w[atom_ctr] = 1.0;
//            atom_ctr++;
//        }
//    }

    for (int i = 0; i < N; ++i) {
        fscanf(points_fp, "%lf %lf %lf %lf", &(sources->x[i]), &(sources->y[i]), &(sources->z[i]), &(sources->q[i]));
        sources->w[i] = 1.0;
    }
    fclose(points_fp);

    
    /* Setting up targets */
    
    targets = malloc(sizeof(struct Particles));
    targets->num = grid_dim[0] * grid_dim[1] * grid_dim[2];
    
    targets->xmin = xyz_limits[0];
    targets->ymin = xyz_limits[2];
    targets->zmin = xyz_limits[4];

    targets->xmax = xyz_limits[1];
    targets->ymax = xyz_limits[3];
    targets->zmax = xyz_limits[5];

    targets->xdim = grid_dim[0];
    targets->ydim = grid_dim[1];
    targets->zdim = grid_dim[2];

    targets->xdd = grid_dd[0];
    targets->ydd = grid_dd[1];
    targets->zdd = grid_dd[2];

/*
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->x));
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->y));
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->z));
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->q));

    for (int k = 0; k < grid_dim[2]; ++k) {
        double zz = xyz_limits[4] + k * grid_dd[2];
        int grid_pos_k = k * grid_dim[1] * grid_dim[0];

        for (int j = 0; j < grid_dim[1]; ++j) {
            double yy = xyz_limits[2] + j * grid_dd[1];
            int grid_pos_jk = grid_pos_k + j * grid_dim[0];

            for (int i = 0; i < grid_dim[0]; ++i) {
                double xx = xyz_limits[0] + i * grid_dd[0];

                targets->x[grid_pos_jk + i] = xx;
                targets->y[grid_pos_jk + i] = yy;
                targets->z[grid_pos_jk + i] = zz;
                targets->q[grid_pos_jk + i] = 1.0;
            }
        }
    }
*/


    if (rank == 0) printf("[random cube example] Setup has finished.\n");

    /* Initializing direct and treedriver runs */

    //targets_sample = malloc(sizeof(struct Particles));

    potential = malloc(sizeof(double) * targets->num);
    potential_direct = malloc(sizeof(double) * targets->num);

    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));


#ifdef OPENACC_ENABLED
    #pragma acc set device_num(rank) device_type(acc_device_nvidia)
    #pragma acc init device_type(acc_device_nvidia)
#endif

    STOP_TIMER(&time_run[0]);
    MPI_Barrier(MPI_COMM_WORLD);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Running direct comparison
    //~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (run_direct == 1) {

// Currently I haven't set up sampling here. So I have to run the whole thing.

//        targets_sample->num = targets->num / slice;
//        targets_sample->x = malloc(targets_sample->num * sizeof(double));
//        targets_sample->y = malloc(targets_sample->num * sizeof(double));
//        targets_sample->z = malloc(targets_sample->num * sizeof(double));
//        targets_sample->q = malloc(targets_sample->num * sizeof(double));
//
//        for (int i = 0; i < targets_sample->num; i++) {
//            targets_sample->x[i] = targets->x[i*slice];
//            targets_sample->y[i] = targets->y[i*slice];
//            targets_sample->z[i] = targets->z[i*slice];
//            targets_sample->q[i] = targets->q[i*slice];
//        }

        if (rank == 0) printf("[random cube example] Running direct comparison...\n");

        START_TIMER(&time_run[1]);
        directdriver(sources, targets, run_params, potential_direct, time_direct);
        STOP_TIMER(&time_run[1]);

        //free(targets_sample->x);
        //free(targets_sample->y);
        //free(targets_sample->z);
        //free(targets_sample->q);
        //free(targets_sample);

    }

    MPI_Barrier(MPI_COMM_WORLD);
    

    //~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Running treecode
    //~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    if (rank == 0) printf("[random cube example] Running treedriver...\n");

    START_TIMER(&time_run[2]);
    treedriver(sources, targets, run_params, potential, time_tree);
    STOP_TIMER(&time_run[2]);

    
    MPI_Barrier(MPI_COMM_WORLD);
    /* Ending total runtime timer */
    STOP_TIMER(&time_run[3]);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Calculate results
    //~~~~~~~~~~~~~~~~~~~~~~~~~~

    Timing_Calculate(time_run_glob, time_tree_glob, time_direct_glob,
                     time_run, time_tree, time_direct);
    Timing_Print(time_run_glob, time_tree_glob, time_direct_glob, run_direct, run_params);
    
    if (run_direct == 1) {
        Accuracy_Calculate(&potential_engy_glob, &potential_engy_direct_glob,
                           &glob_inf_err, &glob_relinf_err, &glob_n2_err, &glob_reln2_err,
                           potential, potential_direct, targets->num, slice);
        Accuracy_Print(potential_engy_glob, potential_engy_direct_glob,
                           glob_inf_err, glob_relinf_err, glob_n2_err, glob_reln2_err, slice);
    }
    
    CSV_Print(N, targets->num, run_params, time_run_glob, time_tree_glob, time_direct_glob,
              potential_engy_glob, potential_engy_direct_glob,
              glob_inf_err, glob_relinf_err, glob_n2_err, glob_reln2_err);


    //~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Cleanup
    //~~~~~~~~~~~~~~~~~~~~~~~~~~

    MPI_Free_mem(sources->x);
    MPI_Free_mem(sources->y);
    MPI_Free_mem(sources->z);
    MPI_Free_mem(sources->q);
    MPI_Free_mem(sources->w);
    free(sources);

    MPI_Free_mem(targets->x);
    MPI_Free_mem(targets->y);
    MPI_Free_mem(targets->z);
    MPI_Free_mem(targets->q);
    free(targets);

    free(potential);
    free(potential_direct);

    RunParams_Free(&run_params);

    MPI_Finalize();

    return 0;
}
