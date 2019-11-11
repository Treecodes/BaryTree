#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <zoltan.h>
#include <time.h>
#include <float.h>

#include "../src/treedriver.h"
#include "../src/particles.h"
#include "../src/tools.h"


const unsigned m = 1664525u;
const unsigned c = 1013904223u;

typedef struct{
  int numGlobalPoints;
  int numMyPoints;
  ZOLTAN_ID_PTR myGlobalIDs;
  double *x;
  double *y;
  double *z;
  double *q;
  double *w;
} MESH_DATA;

typedef struct{
  ZOLTAN_ID_TYPE myGlobalID;
  double x;
  double y;
  double z;
  double q;
  double w;
} SINGLE_MESH_DATA;

/* Application defined query functions */

static int get_number_of_objects(void *data, int *ierr);
static void get_object_list(void *data, int sizeGID, int sizeLID,
            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                  int wgt_dim, float *obj_wgts, int *ierr);
static int get_num_geometry(void *data, int *ierr);
static void get_geometry_list(void *data, int sizeGID, int sizeLID,
             int num_obj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
             int num_dim, double *geom_vec, int *ierr);
static void ztn_pack(void *data, int num_gid_entries, int num_lid_entries,
              ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
              int dest, int size, char *buf, int *ierr);
static void ztn_unpack(void *data, int num_gid_entries,
                ZOLTAN_ID_PTR global_id,
                int size, char *buf, int *ierr);
static int ztn_obj_size(void *data, int num_gid_entries, int num_lid_entries,
    ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr);

int main(int argc, char **argv)
{
    //run parameters
    int N, tree_type, order, max_per_leaf, max_per_batch, run_direct_comparison;
    double kappa, theta; 
    char *kernelName = NULL;
    char *singularityHandling = NULL;
    char *approximationName = NULL;

    N = atoi(argv[1]);
    order = atoi(argv[2]);
    theta = atof(argv[3]);
    max_per_leaf = atoi(argv[4]);
    max_per_batch = atoi(argv[5]);
    kernelName = argv[6];
    kappa = atof(argv[7]);
    singularityHandling = argv[8];
    approximationName = argv[9];
    tree_type = atoi(argv[10]);
    run_direct_comparison = atoi(argv[11]);

//    printf("Read in command line args...\n");
//    mpirun -n ${NP} zoltan_example_cpu $N $ORDER $THETA $CLUSTERSIZE $BATCHSIZE $KERNEL $KAPPA $TREETYPE $COMPAREDIRECT

    int rc, rank, numProcs;
//    printf("Initializing MPI.\n");
    MPI_Init(&argc, &argv);
//    printf("Initialized MPI.\n");
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//    printf("Set rank.\n");
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
//    printf("Set numProcs.\n");

    if (rank == 0) fprintf(stderr,"Beginning BaryTree with %d ranks.\n", numProcs);


    float ver;
    struct Zoltan_Struct *zz;
    int changes, numGidEntries, numLidEntries, numImport, numExport;
    ZOLTAN_ID_PTR importGlobalGids, importLocalGids, exportGlobalGids, exportLocalGids; 
    int *importProcs, *importToPart, *exportProcs, *exportToPart;
    int *parts;
    MESH_DATA mySources;


    struct particles *sources = NULL;
    struct particles *targets = NULL;
    int *particleOrder = NULL;
    double *potential = NULL, *potential_direct = NULL;
    double potential_engy = 0, potential_engy_glob = 0;
    double potential_engy_direct = 0, potential_engy_direct_glob = 0;

    /* variables for date-time calculation */
    double time_run[3], time_tree[9], time_direct;
    double time_run_glob[3][3], time_tree_glob[3][9], time_direct_glob[3];
    double time1, time2;


//    printf("Initializing Zoltan.\n");
    if (Zoltan_Initialize(argc, argv, &ver) != ZOLTAN_OK) {
        printf("Zoltan failed to initialize. Exiting.\n");
        MPI_Finalize();
        exit(0);
    }
//    printf("Initialized Zoltan.\n");

    time1 = MPI_Wtime();

    mySources.numGlobalPoints = N * numProcs;
    mySources.numMyPoints = N;
    mySources.x = malloc(N*sizeof(double));
    mySources.y = malloc(N*sizeof(double));
    mySources.z = malloc(N*sizeof(double));
    mySources.q = malloc(N*sizeof(double));
    mySources.w = malloc(N*sizeof(double));
    mySources.myGlobalIDs = (ZOLTAN_ID_TYPE *)malloc(sizeof(ZOLTAN_ID_TYPE) * N);

    time_t t = time(NULL);
    unsigned t_hashed = (unsigned) t;
    t_hashed = m*t_hashed + c;
//    srand(t_hashed ^ rank);

    for (int i = 0; i < N; ++i) {
        mySources.x[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        mySources.y[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        mySources.z[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        double r = sqrt(mySources.x[i]*mySources.x[i] + mySources.y[i]*mySources.y[i] + mySources.z[i]*mySources.z[i]);
//        mySources.q[i] = exp(-10*r);
        mySources.q[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
//        mySources.q[i] = ((double)rand()/(double)(RAND_MAX));
//        mySources.w[i] = 1.;
//        mySources.w[i] = (i%5)*1.0;
        mySources.w[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;

        mySources.myGlobalIDs[i] = (ZOLTAN_ID_TYPE)(rank*N + i);
    }

    printf("\n\n\nx,y,z [0] = %f,%f,%f\n\n\n",mySources.x[0],mySources.y[0],mySources.z[0] );

    zz = Zoltan_Create(MPI_COMM_WORLD);

    /* General parameters */

    Zoltan_Set_Param(zz, "DEBUG_LEVEL", "0");
    Zoltan_Set_Param(zz, "LB_METHOD", "RCB");
    Zoltan_Set_Param(zz, "NUM_GID_ENTRIES", "1"); 
    Zoltan_Set_Param(zz, "NUM_LID_ENTRIES", "1");
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", "0");
    Zoltan_Set_Param(zz, "RETURN_LISTS", "ALL");
    Zoltan_Set_Param(zz, "AUTO_MIGRATE", "TRUE"); 

    /* RCB parameters */

    Zoltan_Set_Param(zz, "RCB_OUTPUT_LEVEL", "0");
    Zoltan_Set_Param(zz, "RCB_RECTILINEAR_BLOCKS", "1"); 

    /* Query functions, to provide geometry to Zoltan */

    Zoltan_Set_Num_Obj_Fn(zz, get_number_of_objects, &mySources);
    Zoltan_Set_Obj_List_Fn(zz, get_object_list, &mySources);
    Zoltan_Set_Num_Geom_Fn(zz, get_num_geometry, &mySources);
    Zoltan_Set_Geom_Multi_Fn(zz, get_geometry_list, &mySources);
    Zoltan_Set_Obj_Size_Fn(zz, ztn_obj_size, &mySources);
    Zoltan_Set_Pack_Obj_Fn(zz, ztn_pack, &mySources);
    Zoltan_Set_Unpack_Obj_Fn(zz, ztn_unpack, &mySources);

    double x_min = minval(mySources.x, mySources.numMyPoints);
    double x_max = maxval(mySources.x, mySources.numMyPoints);
    double y_min = minval(mySources.y, mySources.numMyPoints);
    double y_max = maxval(mySources.y, mySources.numMyPoints);
    double z_min = minval(mySources.z, mySources.numMyPoints);
    double z_max = maxval(mySources.z, mySources.numMyPoints);

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

    x_min = minval(mySources.x, mySources.numMyPoints);
    x_max = maxval(mySources.x, mySources.numMyPoints);
    y_min = minval(mySources.y, mySources.numMyPoints);
    y_max = maxval(mySources.y, mySources.numMyPoints);
    z_min = minval(mySources.z, mySources.numMyPoints);
    z_max = maxval(mySources.z, mySources.numMyPoints);

    if (rc != ZOLTAN_OK) {
        printf("Error! Zoltan has failed. Exiting. \n");
        MPI_Finalize();
        Zoltan_Destroy(&zz);
        exit(0);
    }

    if (rank == 0) fprintf(stderr,"Zoltan load balancing has finished.\n");

    sources = malloc(sizeof(struct particles));
    targets = malloc(sizeof(struct particles));
    potential = malloc(sizeof(double) * mySources.numMyPoints);
    potential_direct = malloc(sizeof(double) * mySources.numMyPoints);
    particleOrder = malloc(sizeof(int) * mySources.numMyPoints);
    for (int i = 0; i < mySources.numMyPoints; i++) particleOrder[i] = i;

    sources->num = mySources.numMyPoints;
    targets->num = mySources.numMyPoints;

    sources->x = mySources.x;
    sources->y = mySources.y;
    sources->z = mySources.z;
    sources->q = mySources.q;
    sources->w = mySources.w;
    sources->order = particleOrder;
 
//    targets->x = mySources.x;
//    targets->y = mySources.y;
//    targets->z = mySources.z;
//    targets->q = mySources.q;
//    targets->order = particleOrder;
    
    MPI_Barrier(MPI_COMM_WORLD);
    targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));
    targets->order = malloc(targets->num*sizeof(double));
    printf("Beginning memcpy()\n");
    memcpy(targets->x, mySources.x, sizeof(mySources.x));
    memcpy(targets->y, mySources.y, sizeof(mySources.y));
    memcpy(targets->z, mySources.z, sizeof(mySources.z));
    memcpy(targets->q, mySources.q, sizeof(mySources.q));
    memcpy(targets->order, particleOrder, sizeof(particleOrder));
    printf("Ending memcpy()\n");

    double *target_quadrature_weights;
    target_quadrature_weights = malloc(sizeof(double) * targets->num);
    for (int i = 0; i < mySources.numMyPoints; i++) target_quadrature_weights[i] = mySources.w[i];

    // Set up kernel
//    printf("Kernel Name: %s\n",kernelName);
    if       (strcmp(singularityHandling,"skipping")==0){
        printf("Zeroing out initial potential.\n");
        for (int i=0; i<targets->num; i++){
            potential[i]=0.0;
            potential_direct[i]=0.0;
        }

    }else if (strcmp(singularityHandling,"subtraction")==0){
        if (strcmp(kernelName,"coulomb")==0){
            for (int i=0; i<targets->num; i++){
                potential[i]=2.0*M_PI*kappa*kappa*targets->q[i];
                potential_direct[i]=2.0*M_PI*kappa*kappa*targets->q[i];
            }
        }else if (strcmp(kernelName,"yukawa")==0){
            for (int i=0; i<targets->num; i++){
                potential[i]=4.0*M_PI*targets->q[i]/kappa/kappa;  // 4*pi*f_t/k**2
                potential_direct[i]=4.0*M_PI*targets->q[i]/kappa/kappa;  // 4*pi*f_t/k**2
            }
        }else{
            if (rank==0) printf("singularityHandling = %s.\n", singularityHandling);
            printf("Invalid option for singularityHandling. Aborting.\n");
            exit(1);
        }



    }else{
        if (rank==0) printf("kernelName = %s.\n", kernelName);
        if (rank==0) printf("Invalid command line argument for kernelName or ... aborting.\n");
        exit(1);
    }



#ifdef OPENACC_ENABLED
    #pragma acc set device_num(rank) device_type(acc_device_nvidia)
    #pragma acc init device_type(acc_device_nvidia)
#endif

    time_run[0] = MPI_Wtime() - time1;

    /* Calling main treecode subroutine to calculate approximate energy */

    MPI_Barrier(MPI_COMM_WORLD);

    if (run_direct_comparison == 1) {
        if (rank == 0) fprintf(stderr,"Running direct comparison...\n");
        time1 = MPI_Wtime();
        treedriver(sources, targets, 0, 0.0, max_per_leaf, max_per_batch,
                   kernelName, kappa, singularityHandling, approximationName, tree_type, potential_direct, &potential_engy_direct, time_tree);
        time_direct = MPI_Wtime() - time1;
    }
    
    if (rank == 0) fprintf(stderr,"Running treedriver...\n");
    time1 = MPI_Wtime();
    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               kernelName, kappa, singularityHandling, approximationName, tree_type, potential, &potential_engy, time_tree);

    if (rank == 0) fprintf(stderr,"Treedriver has finished.\n");
    
    time_run[1] = MPI_Wtime() - time1;
    time_run[2] = time_run[0] + time_run[1];
    
    MPI_Barrier(MPI_COMM_WORLD);

    
    /* Reducing values to root process */
    MPI_Reduce(time_tree, &time_tree_glob[0], 9, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[1], 9, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[2], 9, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(time_run, &time_run_glob[0], 3, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_run, &time_run_glob[1], 3, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_run, &time_run_glob[2], 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&time_direct, &time_direct_glob[0], 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_glob[1], 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time_direct, &time_direct_glob[2], 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(&potential_engy, &potential_engy_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&potential_engy_direct, &potential_engy_direct_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    
    if (rank == 0)
    {
    
        double min_percent = 100. / time_run_glob[0][2];
        double max_percent = 100. / time_run_glob[1][2];
        double avg_percent = 100. / time_run_glob[2][2];

        /* Printing direct and treecode time calculations: */
        printf("\n\nTreecode timing summary (all times in seconds)...\n\n");
        printf("                                       Avg                           Min                         Max\n");
        printf("|    Total time......................  %9.3e s    (100.00%)      %9.3e s    (100.00%)    %9.3e s    (100.00%) \n",
                     time_run_glob[2][2]/numProcs, time_run_glob[0][2], time_run_glob[1][2]);
        printf("|    |\n");
        printf("|    |....Pre-process................  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n",
                     time_run_glob[2][0]/numProcs, time_run_glob[2][0] * avg_percent,
                     time_run_glob[0][0], time_run_glob[0][0] * min_percent,
                     time_run_glob[1][0], time_run_glob[1][0] * max_percent);
        printf("|    |....Treedriver.................  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n",
                     time_run_glob[2][1]/numProcs, time_run_glob[2][1] * avg_percent,
                     time_run_glob[0][1], time_run_glob[0][1] * min_percent,
                     time_run_glob[1][1], time_run_glob[1][1] * max_percent);
        printf("|         |\n");
        printf("|         |....Setup.................  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n",
                     time_tree_glob[2][0]/numProcs, time_tree_glob[2][0] * avg_percent,
                     time_tree_glob[0][0], time_tree_glob[0][0] * min_percent,
                     time_tree_glob[1][0], time_tree_glob[1][0] * max_percent);
        printf("|              |\n");
        printf("|              |....Build tree.......  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n",
                     time_tree_glob[2][1]/numProcs, time_tree_glob[2][1] * avg_percent,
                     time_tree_glob[0][1], time_tree_glob[0][1] * min_percent,
                     time_tree_glob[1][1], time_tree_glob[1][1] * max_percent);
        printf("|              |....Build array......  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n",
                     time_tree_glob[2][2]/numProcs, time_tree_glob[2][2] * avg_percent,
                     time_tree_glob[0][2], time_tree_glob[0][2] * min_percent,
                     time_tree_glob[1][2], time_tree_glob[1][2] * max_percent);
        printf("|              |....Build batches....  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n",
                     time_tree_glob[2][3]/numProcs, time_tree_glob[2][3] * avg_percent,
                     time_tree_glob[0][3], time_tree_glob[0][3] * min_percent,
                     time_tree_glob[1][3], time_tree_glob[1][3] * max_percent);
        printf("|              |....Fill clusters....  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n",
                     time_tree_glob[2][4]/numProcs, time_tree_glob[2][4] * avg_percent,
                     time_tree_glob[0][4], time_tree_glob[0][4] * min_percent,
                     time_tree_glob[1][4], time_tree_glob[1][4] * max_percent);
        printf("|              |....Build LET........  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n",
                     time_tree_glob[2][5]/numProcs, time_tree_glob[2][5] * avg_percent,
                     time_tree_glob[0][5], time_tree_glob[0][5] * min_percent,
                     time_tree_glob[1][5], time_tree_glob[1][5] * max_percent);
        printf("|              |....Build lists......  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n",
                     time_tree_glob[2][6]/numProcs, time_tree_glob[2][6] * avg_percent,
                     time_tree_glob[0][6], time_tree_glob[0][6] * min_percent,
                     time_tree_glob[1][6], time_tree_glob[1][6] * max_percent);
        printf("|         |\n");
        printf("|         |....Computation...........  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n",
                     time_tree_glob[2][7]/numProcs, time_tree_glob[2][7] * avg_percent,
                     time_tree_glob[0][7], time_tree_glob[0][7] * min_percent,
                     time_tree_glob[1][7], time_tree_glob[1][7] * max_percent);
        printf("|         |....Cleanup...............  %9.3e s    (%6.2f%)      %9.3e s    (%6.2f%)    %9.3e s    (%6.2f%) \n\n",
                     time_tree_glob[2][8]/numProcs, time_tree_glob[2][8] * avg_percent,
                     time_tree_glob[0][8], time_tree_glob[0][8] * min_percent,
                     time_tree_glob[1][8], time_tree_glob[1][8] * max_percent);
        
        printf("               Tree potential energy:  %f\n\n", potential_engy_glob);

        if (run_direct_comparison == 1) {
            printf("\nDirect timing summary (all times in seconds)...\n\n");
            printf("                                       Avg                           Min                         Max\n");
            printf("|    Total time......................  %9.3e s    (100.00%)      %9.3e s    (100.00%)    %9.3e s    (100.00%) \n\n",
                         time_direct_glob[2]/numProcs, time_direct_glob[0], time_direct_glob[1]);
            printf("             Direct potential energy:  %f\n\n\n", potential_engy_direct_glob);
            printf("  Absolute error for total potential:  %e\n",
                   fabs(potential_engy_glob-potential_engy_direct_glob));
            printf("  Relative error for total potential:  %e\n\n",
                   fabs((potential_engy_glob-potential_engy_direct_glob)/potential_engy_direct_glob));
        }
    }

    if (run_direct_comparison == 1) {
        double glob_reln2_err, glob_relinf_err, glob_n2_err, glob_inf_err;
        double inferr = 0.0, relinferr = 0.0, n2err = 0.0, reln2err = 0.0;
        double temp;

        for (int j = 0; j < targets->num; ++j) {
            temp = fabs(potential_direct[j] - potential[j]);

            if (temp >= inferr) inferr = temp;

            if (fabs(potential_direct[j]) >= relinferr)
                relinferr = fabs(potential_direct[j]);

            n2err = n2err + pow(potential_direct[j] - potential[j], 2.0) * target_quadrature_weights[j];
            reln2err = reln2err + pow(potential_direct[j], 2.0) * target_quadrature_weights[j];
        }

        MPI_Reduce(&reln2err, &glob_reln2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&n2err, &glob_n2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&relinferr, &glob_relinf_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        MPI_Reduce(&inferr, &glob_inf_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    
        if (rank == 0) {
            glob_reln2_err = sqrt(glob_n2_err / glob_reln2_err);
            glob_n2_err = sqrt(glob_n2_err);
            glob_relinf_err = glob_inf_err / glob_relinf_err;
            printf("Relative inf norm error in potential:  %e \n", glob_relinf_err);
            printf("  Relative 2 norm error in potential:  %e \n\n", glob_reln2_err);
        }
    }

    //for (int i = 0; i < mySources.numMyPoints; i++) {
    //    fprintf(stderr,"Rank %d, particle %d: %f, %f, %f, %d\n", rank, i, mySources.x[i], mySources.y[i], mySources.z[i], mySources.myGlobalIDs[i]);
    //}


    /******************************************************************
     ** Free the arrays allocated by Zoltan_LB_Partition, and free
     ** the storage allocated for the Zoltan structure.
     ******************************************************************/

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
    free(mySources.myGlobalIDs);
    free(particleOrder);
    free(potential);
    free(potential_direct);
    free(sources);
    free(targets);

    MPI_Finalize();

    return 0;
}


static int get_number_of_objects(void *data, int *ierr)
{
    MESH_DATA *mesh = (MESH_DATA *)data;
    *ierr = ZOLTAN_OK;
    return mesh->numMyPoints;
}

static void get_object_list(void *data, int sizeGID, int sizeLID,
                            ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                            int wgt_dim, float *obj_wgts, int *ierr)
{
    int i;
    MESH_DATA *mesh = (MESH_DATA *)data;
    *ierr = ZOLTAN_OK;

  /* In this example, return the IDs of our objects, but no weights.
 *    * Zoltan will assume equally weighted objects.
 *       */

    for (i = 0; i < mesh->numMyPoints; i++) {
        globalID[i] = mesh->myGlobalIDs[i];
        localID[i] = i;
    }
}

static int get_num_geometry(void *data, int *ierr)
{
    *ierr = ZOLTAN_OK;
    return 3;
}

static void get_geometry_list(void *data, int sizeGID, int sizeLID, int num_obj,
                              ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
                              int num_dim, double *geom_vec, int *ierr)
{
    int i;

    MESH_DATA *mesh = (MESH_DATA *)data;

    if ( (sizeGID != 1) || (sizeLID != 1) || (num_dim != 3)) {
        *ierr = ZOLTAN_FATAL;
        return;
    }

    *ierr = ZOLTAN_OK;

    for (i = 0;  i < num_obj ; i++){
        geom_vec[3*i] = (double)mesh->x[i];
        geom_vec[3*i + 1] = (double)mesh->y[i];
        geom_vec[3*i + 2] = (double)mesh->z[i];
    } 

    return;
}

static void ztn_pack(void *data, int num_gid_entries, int num_lid_entries,
              ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id,
              int dest, int size, char *buf, int *ierr) {

    SINGLE_MESH_DATA *mesh_single = (SINGLE_MESH_DATA *)buf;
    MESH_DATA *mesh = (MESH_DATA *)data;

    mesh_single->x = mesh->x[(*local_id)];
    mesh_single->y = mesh->y[(*local_id)];
    mesh_single->z = mesh->z[(*local_id)];
    mesh_single->q = mesh->q[(*local_id)];
    mesh_single->w = mesh->w[(*local_id)];
    mesh_single->myGlobalID = mesh->myGlobalIDs[(*local_id)];

    mesh->myGlobalIDs[(*local_id)] = (ZOLTAN_ID_TYPE)(-1); // Mark local particle as exported

    return;
}

static void ztn_unpack(void *data, int num_gid_entries,
                ZOLTAN_ID_PTR global_id,
                int size, char *buf, int *ierr) {

    SINGLE_MESH_DATA *mesh_single = (SINGLE_MESH_DATA *)buf;
    MESH_DATA *mesh = (MESH_DATA *)data;

    mesh->numMyPoints += 1;

    mesh->myGlobalIDs = (ZOLTAN_ID_TYPE *)realloc(mesh->myGlobalIDs,
                        sizeof(ZOLTAN_ID_TYPE) * mesh->numMyPoints);
    mesh->x = (double *)realloc(mesh->x, sizeof(double) * mesh->numMyPoints);
    mesh->y = (double *)realloc(mesh->y, sizeof(double) * mesh->numMyPoints);
    mesh->z = (double *)realloc(mesh->z, sizeof(double) * mesh->numMyPoints);
    mesh->q = (double *)realloc(mesh->q, sizeof(double) * mesh->numMyPoints);
    mesh->w = (double *)realloc(mesh->w, sizeof(double) * mesh->numMyPoints);

    mesh->x[mesh->numMyPoints-1] = mesh_single->x;
    mesh->y[mesh->numMyPoints-1] = mesh_single->y;
    mesh->z[mesh->numMyPoints-1] = mesh_single->z;
    mesh->q[mesh->numMyPoints-1] = mesh_single->q;
    mesh->w[mesh->numMyPoints-1] = mesh_single->w;
    mesh->myGlobalIDs[mesh->numMyPoints-1] = mesh_single->myGlobalID;

    return;
}

static int ztn_obj_size(void *data, int num_gid_entries, int num_lid_entries, 
    ZOLTAN_ID_PTR global_id, ZOLTAN_ID_PTR local_id, int *ierr)
{
    return sizeof(SINGLE_MESH_DATA);
}
