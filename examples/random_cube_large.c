#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <mpi.h>

#include "../src/utilities/tools.h"

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
    int N, M, run_direct, slice;
    struct RunParams *run_params = NULL;
    int sample_size = 10000;

    FILE *fp = fopen(argv[1], "r");

    Parse_Params(fp, &run_params, &N, &M, &run_direct, &slice);


    int rc, rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    double timebeg = MPI_Wtime();

    if (rank == 0) fprintf(stderr,"Beginning BaryTree with %d ranks.\n", numProcs);


    float ver;
    struct Zoltan_Struct *zz;
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids; 
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    int *parts;
    MESH_DATA mySources;


    struct Particles *sources = NULL;
    struct Particles *targets = NULL;
    struct Particles *targets_sample = NULL;
    double *potential = NULL, *potential_direct = NULL;
    double potential_engy = 0, potential_engy_glob = 0;
    double potential_engy_direct = 0, potential_engy_direct_glob = 0;

    /* variables for date-time calculation */
    double time_run[4], time_tree[13], time_direct[4];
    double time_run_glob[3][4], time_tree_glob[3][13], time_direct_glob[3][4];
    double time1, time2;


    if (Zoltan_Initialize(argc, argv, &ver) != ZOLTAN_OK) {
        if (rank == 0) printf("Zoltan failed to initialize. Exiting.\n");
        MPI_Finalize();
        exit(0);
    }


    time1 = MPI_Wtime();

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
    srand(t_hashed ^ rank);
    srand(1);

    for (int i = 0; i < sample_size; ++i) {
        mySources.x[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        mySources.y[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        mySources.z[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        mySources.q[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        mySources.w[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
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
            mySources.q[i] = mySources.q[mySources.numMyPoints-1];
            mySources.w[i] = mySources.w[mySources.numMyPoints-1];
            mySources.myGlobalIDs[i] = mySources.myGlobalIDs[mySources.numMyPoints-1];
            mySources.numMyPoints--; 
        } else {
          i++;
        }
    }

    if (rc != ZOLTAN_OK) {
        printf("Error! Zoltan has failed. Exiting. \n");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

    double xmin = minval(mySources.x, mySources.numMyPoints);
    double ymin = minval(mySources.y, mySources.numMyPoints);
    double zmin = minval(mySources.z, mySources.numMyPoints);
    double xmax = maxval(mySources.x, mySources.numMyPoints);
    double ymax = maxval(mySources.y, mySources.numMyPoints);
    double zmax = maxval(mySources.z, mySources.numMyPoints);


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

    if (rank == 0) fprintf(stderr,"Zoltan load balancing has finished.\n");


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


    //MPI-allocated source and target arrays for RMA use
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->x));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->y));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->z));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->q));
    MPI_Alloc_mem(sources->num * sizeof(double), MPI_INFO_NULL, &(sources->w));

    //Generating sources and targets based on Zoltan bounding box
    for (int i = 0; i < sources->num; ++i) {
        sources->x[i] = ((double)rand()/(double)(RAND_MAX)) * (xmax-xmin) + xmin;
        sources->y[i] = ((double)rand()/(double)(RAND_MAX)) * (ymax-ymin) + ymin;
        sources->z[i] = ((double)rand()/(double)(RAND_MAX)) * (zmax-zmin) + zmin;
        sources->q[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->w[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
    }

    //MPI-allocated target arrays for RMA use
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->x));
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->y));
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->z));
    MPI_Alloc_mem(targets->num * sizeof(double), MPI_INFO_NULL, &(targets->q));

    //Generating targets based on Zoltan bounding box
    for (int i = 0; i < targets->num; ++i) {
        targets->x[i] = ((double)rand()/(double)(RAND_MAX)) * (xmax-xmin) + xmin;
        targets->y[i] = ((double)rand()/(double)(RAND_MAX)) * (ymax-ymin) + ymin;
        targets->z[i] = ((double)rand()/(double)(RAND_MAX)) * (zmax-zmin) + zmin;
        targets->q[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
    }

#ifdef OPENACC_ENABLED
    #pragma acc set device_num(rank) device_type(acc_device_nvidia)
    #pragma acc init device_type(acc_device_nvidia)
#endif

    time_run[0] = MPI_Wtime() - time1;

    /* Calling main treecode subroutine to calculate approximate energy */

    MPI_Barrier(MPI_COMM_WORLD);

    if (run_direct == 1) {

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

        if (rank == 0) fprintf(stderr,"Running direct comparison...\n");

        time1 = MPI_Wtime();
        directdriver(sources, targets_sample, run_params, potential_direct, time_direct);
        time_run[1] = MPI_Wtime() - time1;

        potential_engy_direct = sum(potential_direct, targets_sample->num);

        free(targets_sample->x);
        free(targets_sample->y);
        free(targets_sample->z);
        free(targets_sample->q);
        free(targets_sample);

    }


    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) fprintf(stderr,"Running treedriver...\n");

    time1 = MPI_Wtime();
    treedriver(sources, targets, run_params, potential, time_tree);
    time_run[2] = MPI_Wtime() - time1;

    potential_engy = sum(potential, targets->num);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    time_run[3] = MPI_Wtime() - timebeg;


    /* Reducing values to root process */
    MPI_Reduce(time_tree, &time_tree_glob[0], 13, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[1], 13, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[2], 13, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(time_run, &time_run_glob[0], 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_run, &time_run_glob[1], 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_run, &time_run_glob[2], 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&time_direct, &time_direct_glob[0], 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_glob[1], 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_glob[2], 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&potential_engy, &potential_engy_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&potential_engy_direct, &potential_engy_direct_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    
    if (rank == 0)
    {
    
        double min_percent = 100. / time_run_glob[0][3];
        double max_percent = 100. / time_run_glob[1][3];
        double avg_percent = 100. / time_run_glob[2][3];

        double min_percent_direct = 100. / time_run_glob[0][1];
        double max_percent_direct = 100. / time_run_glob[1][1];
        double avg_percent_direct = 100. / time_run_glob[2][1];

        double min_percent_tree = 100. / time_run_glob[0][2];
        double max_percent_tree = 100. / time_run_glob[1][2];
        double avg_percent_tree = 100. / time_run_glob[2][2];

        /* Printing direct and treecode time calculations: */
        printf("\n\nTreecode timing summary (all times in seconds)...\n\n");
        printf("                                       Max                           Avg                          Max/Min\n");
        printf("|    Total time......................  %9.3e s    (100.00%%)      %9.3e s    (100.00%%)    %8.3f \n",
                     time_run_glob[1][3], time_run_glob[2][3]/numProcs, time_run_glob[1][3]/time_run_glob[0][3]);
        printf("|    |\n");
        printf("|    |....Pre-process................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_run_glob[1][0],          time_run_glob[1][0] * max_percent,
                     time_run_glob[2][0]/numProcs, time_run_glob[2][0] * avg_percent,
                     time_run_glob[1][0]/time_run_glob[0][0]);

        if (run_direct == 1) {
        printf("|    |....Directdriver...............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_run_glob[1][1],          time_run_glob[1][1] * max_percent,
                     time_run_glob[2][1]/numProcs, time_run_glob[2][1] * avg_percent,
                     time_run_glob[1][1]/time_run_glob[0][1]);
        }
        printf("|    |....Treedriver.................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n\n\n",
                     time_run_glob[1][2],          time_run_glob[1][2] * max_percent,
                     time_run_glob[2][2]/numProcs, time_run_glob[2][2] * avg_percent,
                     time_run_glob[1][2]/time_run_glob[0][2]);


        if (run_direct == 1) {
        printf("|    Directdriver....................  %9.3e s    (100.00%%)      %9.3e s    (100.00%%)    %8.3f \n",
                     time_run_glob[1][1], time_run_glob[2][1]/numProcs, time_run_glob[1][1]/time_run_glob[0][1]);

        printf("|    |\n");

        if (numProcs > 1) {
        printf("|    |....Communicate................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_direct_glob[1][0],          time_direct_glob[1][0] * max_percent_direct,
                     time_direct_glob[2][0]/numProcs, time_direct_glob[2][0] * avg_percent_direct,
                     time_direct_glob[1][0]/time_direct_glob[0][0]);
        }

        printf("|    |....Compute local..............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_direct_glob[1][2],          time_direct_glob[1][2] * max_percent_direct,
                     time_direct_glob[2][2]/numProcs, time_direct_glob[2][2] * avg_percent_direct,
                     time_direct_glob[1][2]/time_direct_glob[0][2]);

        if (numProcs > 1) {
        printf("|    |....Compute remote.............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_direct_glob[1][1],          time_direct_glob[1][1] * max_percent_direct,
                     time_direct_glob[2][1]/numProcs, time_direct_glob[2][1] * avg_percent_direct,
                     time_direct_glob[1][1]/time_direct_glob[0][1]);
        }

        if (run_params->singularity == SUBTRACTION) {
        printf("|    |....Correct potential..........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n\n\n",
                     time_direct_glob[1][3],          time_direct_glob[1][3] * max_percent_direct,
                     time_direct_glob[2][3]/numProcs, time_direct_glob[2][3] * avg_percent_direct,
                     time_direct_glob[1][3]/time_direct_glob[0][3]);
        } else {
        printf("\n\n");
        }
        }


        printf("|    Treedriver......................  %9.3e s    (100.00%%)      %9.3e s    (100.00%%)    %8.3f \n",
                     time_run_glob[1][2], time_run_glob[2][2]/numProcs, time_run_glob[1][2]/time_run_glob[0][2]);
        printf("|    |\n");
        printf("|    |....Build local tree...........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_tree_glob[1][0],          time_tree_glob[1][0] * max_percent_tree,
                     time_tree_glob[2][0]/numProcs, time_tree_glob[2][0] * avg_percent_tree,
                     time_tree_glob[1][0]/time_tree_glob[0][0]);
        printf("|    |....Build local batches........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_tree_glob[1][1],          time_tree_glob[1][1] * max_percent_tree,
                     time_tree_glob[2][1]/numProcs, time_tree_glob[2][1] * avg_percent_tree,
                     time_tree_glob[1][1]/time_tree_glob[0][1]);
        printf("|    |....Build local clusters.......  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_tree_glob[1][2],          time_tree_glob[1][2] * max_percent_tree,
                     time_tree_glob[2][2]/numProcs, time_tree_glob[2][2] * avg_percent_tree,
                     time_tree_glob[1][2]/time_tree_glob[0][2]);

        if (numProcs > 1) {
        printf("|    |....Build LET..................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_tree_glob[1][3],          time_tree_glob[1][3] * max_percent_tree,
                     time_tree_glob[2][3]/numProcs, time_tree_glob[2][3] * avg_percent_tree,
                     time_tree_glob[1][3]/time_tree_glob[0][3]);
        }

        printf("|    |....Build local lists..........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_tree_glob[1][4],          time_tree_glob[1][4] * max_percent_tree,
                     time_tree_glob[2][4]/numProcs, time_tree_glob[2][4] * avg_percent_tree,
                     time_tree_glob[1][4]/time_tree_glob[0][4]);
        printf("|    |....Compute local..............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_tree_glob[1][5],          time_tree_glob[1][5] * max_percent_tree,
                     time_tree_glob[2][5]/numProcs, time_tree_glob[2][5] * avg_percent_tree,
                     time_tree_glob[1][5]/time_tree_glob[0][5]);

        if (numProcs > 1) {
        printf("|    |....Build remote lists.........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_tree_glob[1][6],          time_tree_glob[1][6] * max_percent_tree,
                     time_tree_glob[2][6]/numProcs, time_tree_glob[2][6] * avg_percent_tree,
                     time_tree_glob[1][6]/time_tree_glob[0][6]);
        printf("|    |....Compute remote.............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_tree_glob[1][7],          time_tree_glob[1][7] * max_percent_tree,
                     time_tree_glob[2][7]/numProcs, time_tree_glob[2][7] * avg_percent_tree,
                     time_tree_glob[1][7]/time_tree_glob[0][7]);
        }

        if (run_params->compute_type == CLUSTER_PARTICLE || run_params->compute_type == CLUSTER_CLUSTER) {
        printf("|    |....Compute cp2................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_tree_glob[1][8],          time_tree_glob[1][8] * max_percent_tree,
                     time_tree_glob[2][8]/numProcs, time_tree_glob[2][8] * avg_percent_tree,
                     time_tree_glob[1][8]/time_tree_glob[0][8]);
        }

        printf("|    |....Correct potential..........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                     time_tree_glob[1][9],          time_tree_glob[1][9] * max_percent_tree,
                     time_tree_glob[2][9]/numProcs, time_tree_glob[2][9] * avg_percent_tree,
                     time_tree_glob[1][9]/time_tree_glob[0][9]);
        printf("|    |....Cleanup....................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n\n",
                     time_tree_glob[1][10],          time_tree_glob[1][10] * max_percent_tree,
                     time_tree_glob[2][10]/numProcs, time_tree_glob[2][10] * avg_percent_tree,
                     time_tree_glob[1][10]/time_tree_glob[0][10]);
        
        if (numProcs > 1) {
        printf("((   |....Total setup................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f ))\n",
                     time_tree_glob[1][11],          time_tree_glob[1][11] * max_percent_tree,
                     time_tree_glob[2][11]/numProcs, time_tree_glob[2][11] * avg_percent_tree,
                     time_tree_glob[1][11]/time_tree_glob[0][11]);
        printf("((   |....Build local clusters.......  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f ))\n",
                     time_tree_glob[1][02],          time_tree_glob[1][02] * max_percent_tree,
                     time_tree_glob[2][02]/numProcs, time_tree_glob[2][02] * avg_percent_tree,
                     time_tree_glob[1][02]/time_tree_glob[0][02]);
        printf("((   |....Total compute..............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f ))\n\n\n",
                     time_tree_glob[1][12],          time_tree_glob[1][12] * max_percent_tree,
                     time_tree_glob[2][12]/numProcs, time_tree_glob[2][12] * avg_percent_tree,
                     time_tree_glob[1][12]/time_tree_glob[0][12]);
        }
        
        printf("               Tree potential energy:  %f\n", potential_engy_glob);

        if (run_direct == 1) {
        printf("             Direct potential energy:  %f\n\n", potential_engy_direct_glob);
        printf("  Absolute error for total potential:  %e\n",
               fabs(potential_engy_glob-potential_engy_direct_glob));
        printf("  Relative error for total potential:  %e\n\n",
               fabs((potential_engy_glob-potential_engy_direct_glob)/potential_engy_direct_glob));
        }
    }

    double glob_reln2_err, glob_relinf_err, glob_n2_err, glob_inf_err;
    if (run_direct == 1) {
        double inferr = 0.0, relinferr = 0.0, n2err = 0.0, reln2err = 0.0;
        double temp;

        for (int j = 0; j < targets->num / slice; j++) {

            temp = fabs(potential_direct[j] - potential[j*slice]);

            if (temp >= inferr) inferr = temp;

            if (fabs(potential_direct[j]) >= relinferr)
                relinferr = fabs(potential_direct[j]);

            n2err = n2err + pow(potential_direct[j] - potential[j*slice], 2.0);
            reln2err = reln2err + pow(potential_direct[j], 2.0);
        }

        MPI_Reduce(&reln2err, &glob_reln2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&n2err, &glob_n2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&relinferr, &glob_relinf_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&inferr, &glob_inf_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
        if (rank == 0) {
            glob_reln2_err = sqrt(fabs(glob_n2_err / glob_reln2_err));
            glob_n2_err = sqrt(fabs(glob_n2_err));
            glob_relinf_err = glob_inf_err / glob_relinf_err;
            printf("Relative inf norm error in potential:  %e \n", glob_relinf_err);
            printf("  Relative 2 norm error in potential:  %e \n\n", glob_reln2_err);
        }
    }


    if (rank == 0) {
        FILE *fp = fopen("out.csv", "a");
        fprintf(fp, "%d,%d,%d,%f,%d,%d,%d,%d,%d,%d,"
                    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                    "%e,%e,%e,%e,%e,%e,%e,%e\n",
            N, M, run_params->interp_order, run_params->theta,
            run_params->max_per_source_leaf, run_params->max_per_target_leaf, run_params->kernel,
            run_params->singularity, run_params->approximation, numProcs, // 1 ends

            time_run_glob[0][0],  time_run_glob[1][0],  // min, max, avg pre-process
            time_run_glob[2][0]/numProcs,

            time_run_glob[0][1],  time_run_glob[1][1],  // min, max, avg directdriver
            time_run_glob[2][1]/numProcs,

            time_direct_glob[0][0], time_direct_glob[1][0],     // min, max, avg communicate
            time_direct_glob[2][0]/numProcs,
            time_direct_glob[0][1], time_direct_glob[1][1],     // min, max, avg compute local
            time_direct_glob[2][1]/numProcs,
            time_direct_glob[0][2], time_direct_glob[1][2],     // min, max, avg compute remote
            time_direct_glob[2][2]/numProcs,
            time_direct_glob[0][3], time_direct_glob[1][3],     // min, max, avg correct potential
            time_direct_glob[2][3]/numProcs, // 2 ends

            time_run_glob[0][2],  time_run_glob[1][2],  // min, max, avg treedriver
            time_run_glob[2][2]/numProcs,

            time_tree_glob[0][0], time_tree_glob[1][0],     // min, max, avg build local tree
            time_tree_glob[2][0]/numProcs,
            time_tree_glob[0][1], time_tree_glob[1][1],     // min, max, avg build local batches
            time_tree_glob[2][1]/numProcs,
            time_tree_glob[0][2], time_tree_glob[1][2],     // min, max, avg fill local clusters
            time_tree_glob[2][2]/numProcs,
            time_tree_glob[0][3], time_tree_glob[1][3],     // min, max, avg build LET
            time_tree_glob[2][3]/numProcs,
            time_tree_glob[0][4], time_tree_glob[1][4],     // min, max, avg build local lists
            time_tree_glob[2][4]/numProcs, // 3 ends

            time_tree_glob[0][5], time_tree_glob[1][5],     // min, max, avg compute local interations 
            time_tree_glob[2][5]/numProcs,
            time_tree_glob[0][6], time_tree_glob[1][6],     // min, max, avg build remote lists
            time_tree_glob[2][6]/numProcs,
            time_tree_glob[0][7], time_tree_glob[1][7],     // min, max, avg compute remote interactions
            time_tree_glob[2][7]/numProcs,
            time_tree_glob[0][8], time_tree_glob[1][8],     // min, max, avg compute CP2 (cluster-particle)
            time_tree_glob[2][8]/numProcs,
            time_tree_glob[0][9], time_tree_glob[1][9],     // min, max, avg correct potential
            time_tree_glob[2][9]/numProcs,
            time_tree_glob[0][10], time_tree_glob[1][10],     // min, max, avg cleanup
            time_tree_glob[2][10]/numProcs,

            time_tree_glob[0][11], time_tree_glob[1][11],   // min, max, avg total setup
            time_tree_glob[2][11]/numProcs,
            time_tree_glob[0][12], time_tree_glob[1][12],   // min, max, avg total cleanup
            time_tree_glob[2][12]/numProcs,

            time_run_glob[0][3],  time_run_glob[1][3],  // min, max, avg total time
            time_run_glob[2][3]/numProcs, // 4 ends

            potential_engy_direct_glob, potential_engy_glob,
            fabs(potential_engy_direct_glob - potential_engy_glob),
            fabs((potential_engy_direct_glob - potential_engy_glob) / potential_engy_direct_glob),
            glob_inf_err, glob_relinf_err, glob_n2_err, glob_reln2_err); // 5 ends
        fclose(fp);
    }


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
