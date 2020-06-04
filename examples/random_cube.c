#include <math.h>
#include <stdlib.h>
#include <stdio.h>
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

#include "zoltan_fns.h"
#include "support_fns.h"

const unsigned mrand = 1664525u;
const unsigned crand = 1013904223u;


int main(int argc, char **argv)
{
    /* MPI initialization */
    int rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    if (rank == 0) printf("[random cube example] Beginning random cube example with %d ranks.\n", numProcs);
    
    /* run parameters */
    int N, M, run_direct, slice;
    double xyz_limits[6];
    DISTRIBUTION distribution;
    PARTITION partition;
    int sample_size = 1000000;
    
    struct RunParams *run_params = NULL;
    
    FILE *fp = fopen(argv[1], "r");
    Params_Parse(fp, &run_params, &N, &M, &run_direct, &slice, xyz_limits, &distribution, &partition);
    
    double xmin = xyz_limits[0], xmax = xyz_limits[1];
    double ymin = xyz_limits[2], ymax = xyz_limits[3];
    double zmin = xyz_limits[4], zmax = xyz_limits[5];

    /* Zoltan variables */
    int rc;
    float ver;
    struct Zoltan_Struct *zz;
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids; 
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    int *parts;
    MESH_DATA mySources;

    /* data structures for BaryTree calculation and comparison */
    struct Particles *sources = NULL;
    struct Particles *targets = NULL;
    struct Particles *targets_sample = NULL;
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
    
    /* Zoltan initialization */
    if (Zoltan_Initialize(argc, argv, &ver) != ZOLTAN_OK) {
        if (rank == 0) printf("[random cube example] Zoltan failed to initialize. Exiting.\n");
        MPI_Finalize();
        exit(0);
    }
    
    mySources.numGlobalPoints = sample_size * numProcs;
    mySources.numMyPoints = sample_size;
    mySources.x = malloc(sample_size * sizeof(double));
    mySources.y = malloc(sample_size * sizeof(double));
    mySources.z = malloc(sample_size * sizeof(double));
    mySources.q = malloc(sample_size * sizeof(double));
    mySources.w = malloc(sample_size * sizeof(double));
    mySources.b = malloc(sample_size * sizeof(double)); // load balancing weights
    mySources.myGlobalIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * sample_size);

    time_t t = time(NULL);
    unsigned t_hashed = (unsigned) t;
    t_hashed = mrand * t_hashed + crand;
    srandom(t_hashed ^ rank);
    //srandom(1);

    for (int i = 0; i < sample_size; ++i) {
        mySources.x[i] = Point_Set_Init(distribution);
        mySources.y[i] = Point_Set_Init(distribution);
        mySources.z[i] = Point_Set_Init(distribution);
        mySources.myGlobalIDs[i] = (ZOLTAN_ID_TYPE)(rank*N + i);

        mySources.b[i] = 1.0; // dummy weighting scheme
    }


    zz = Zoltan_Create(MPI_COMM_WORLD);

    /* General parameters */

    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "1");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
    Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE"); 

    /* RCB parameters */

    Zoltan_Set_Param(zz, "KEEP_CUTS", "1");
    Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
    Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1"); 

    /* Query functions, to provide geometry to Zoltan */

    Zoltan_Set_Num_Obj_Fn(zz, ztn_get_number_of_objects, &mySources);
    Zoltan_Set_Obj_List_Fn(zz, ztn_get_object_list, &mySources);
    Zoltan_Set_Num_Geom_Fn(zz, ztn_get_num_geometry, &mySources);
    Zoltan_Set_Geom_Multi_Fn(zz, ztn_get_geometry_list, &mySources);
    Zoltan_Set_Obj_Size_Fn(zz, ztn_obj_size, &mySources);
    Zoltan_Set_Pack_Obj_Fn(zz, ztn_pack, &mySources);
    Zoltan_Set_Unpack_Obj_Fn(zz, ztn_unpack, &mySources);

    rc = Zoltan_LB_Partition(zz, /* input (all remaining fields are output) */
        &changes,        /* 1 if partitioning was changed, 0 otherwise */ 
        &numGidEntries,  /* Number of integers used for a global ID */
        &numLidEntries,  /* Number of integers used for a local ID */
        &numImport,      /* Number of vertices to be sent to me */
        &importGlobalGids,  /* Global IDs of vertices to be sent to me */
        &importLocalGids,   /* Local IDs of vertices to be sent to me */
        &importProcs,    /* Process rank for source of each incoming vertex */
        &importToPart,   /* New partition for each incoming vertex */
        &numExport,      /* Number of vertices I must send to other processes*/
        &exportGlobalGids,  /* Global IDs of the vertices I must send */
        &exportLocalGids,   /* Local IDs of the vertices I must send */
        &exportProcs,    /* Process to which I send each of the vertices */
        &exportToPart);  /* Partition to which each vertex will belong */

    int i = 0;
    while (i < mySources.numMyPoints) {
        if ((int)mySources.myGlobalIDs[i] < 0) {
            mySources.x[i] = mySources.x[mySources.numMyPoints-1];
            mySources.y[i] = mySources.y[mySources.numMyPoints-1];
            mySources.z[i] = mySources.z[mySources.numMyPoints-1];
            mySources.myGlobalIDs[i] = mySources.myGlobalIDs[mySources.numMyPoints-1];
            mySources.numMyPoints--; 
        } else {
          i++;
        }
    }

    if (rc != ZOLTAN_OK) {
        printf("[random cube example] Error! Zoltan has failed. Exiting. \n");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

    double zz_bound_x_min = minval(mySources.x, mySources.numMyPoints);
    double zz_bound_y_min = minval(mySources.y, mySources.numMyPoints);
    double zz_bound_z_min = minval(mySources.z, mySources.numMyPoints);
    double zz_bound_x_max = maxval(mySources.x, mySources.numMyPoints);
    double zz_bound_y_max = maxval(mySources.y, mySources.numMyPoints);
    double zz_bound_z_max = maxval(mySources.z, mySources.numMyPoints);


    Zoltan_LB_Free_Part(&importGlobalGids, &importLocalGids, 
                        &importProcs, &importToPart);
    Zoltan_LB_Free_Part(&exportGlobalGids, &exportLocalGids, 
                        &exportProcs, &exportToPart);
    Zoltan_Destroy(&zz);

    free(mySources.x);
    free(mySources.y);
    free(mySources.z);
    free(mySources.q);
    free(mySources.w);
    free(mySources.b);
    free(mySources.myGlobalIDs);

    if (rank == 0) printf("[random cube example] Zoltan load balancing has finished.\n");


    sources = malloc(sizeof(struct Particles));
    sources->num = N;

    targets = malloc(sizeof(struct Particles));
    targets->num = M;
    potential = malloc(targets->num * sizeof(double));
    memset(potential, 0, targets->num * sizeof(double));

    targets_sample = malloc(sizeof(struct Particles));
    targets_sample->num = targets->num / slice;
    potential_direct = malloc(targets_sample->num * sizeof(double));
    memset(potential_direct, 0, targets_sample->num * sizeof(double));


    /* MPI-allocated source and target arrays for RMA use */
    
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->x));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->y));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->z));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->q));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->w));

    /* Generating sources and targets based on Zoltan bounding box */
    
    for (int i = 0; i < sources->num; ++i) {
        sources->x[i] = Point_Set(distribution, zz_bound_x_min, zz_bound_x_max) * (xmax-xmin) + xmin;
        sources->y[i] = Point_Set(distribution, zz_bound_y_min, zz_bound_y_max) * (ymax-ymin) + ymin;
        sources->z[i] = Point_Set(distribution, zz_bound_z_min, zz_bound_z_max) * (zmax-zmin) + zmin;

        sources->q[i] = Point_Set(UNIFORM, -1., 1.);
        sources->w[i] = Point_Set(UNIFORM, -1., 1.);
    }

/*
    char points_file[256];
    sprintf(points_file, "points_rank_%d.csv", rank);
    FILE *points_fp = fopen(points_file, "w");
    for (int i = 0; i < sources->num; ++i) {
        fprintf(points_fp, "%e, %e, %e\n", sources->x[i], sources->y[i], sources->z[i]);
    }
    fclose(points_fp);
*/

    /* MPI-allocated target arrays for RMA use */
    
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->x));
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->y));
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->z));
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->q));

    /* Generating targets based on Zoltan bounding box */
    
    for (int i = 0; i < targets->num; ++i) {
        targets->x[i] = Point_Set(distribution, zz_bound_x_min, zz_bound_x_max) * (xmax-xmin) + xmin;
        targets->y[i] = Point_Set(distribution, zz_bound_y_min, zz_bound_y_max) * (ymax-ymin) + ymin;
        targets->z[i] = Point_Set(distribution, zz_bound_z_min, zz_bound_z_max) * (zmax-zmin) + zmin;
        targets->q[i] = Point_Set(UNIFORM, -1., 1.);
    }

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
    
        if (rank == 0) printf("[random cube example] Running direct comparison...\n");

        targets_sample->x = malloc(targets_sample->num * sizeof(double));
        targets_sample->y = malloc(targets_sample->num * sizeof(double));
        targets_sample->z = malloc(targets_sample->num * sizeof(double));
        targets_sample->q = malloc(targets_sample->num * sizeof(double));

        for (int i = 0; i < targets_sample->num; i++) {
            targets_sample->x[i] = targets->x[i*slice];
            targets_sample->y[i] = targets->y[i*slice];
            targets_sample->z[i] = targets->z[i*slice];
            targets_sample->q[i] = targets->q[i*slice];
        }

        START_TIMER(&time_run[1]);
        directdriver(sources, targets_sample, run_params, potential_direct, time_direct);
        STOP_TIMER(&time_run[1]);

        free(targets_sample->x);
        free(targets_sample->y);
        free(targets_sample->z);
        free(targets_sample->q);
        free(targets_sample);

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
    
    CSV_Print(N, M, run_params, time_run_glob, time_tree_glob, time_direct_glob,
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
