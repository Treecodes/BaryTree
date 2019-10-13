#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <zoltan.h>

#include "../src/treedriverWrapper.h"

typedef struct{
  int numGlobalPoints;
  int numMyPoints;
  ZOLTAN_ID_PTR myGlobalIDs;
  double *x;
  double *y;
  double *z;
} MESH_DATA;

int main(int argc, char **argv)
{
    //run parameters
    int N = 100000;
    int M = 100000;

    int pot_type = 0;
    double kappa = 0.;

    int order = 5;
    double theta = 0.5;

    int max_per_leaf = 500;
    int max_per_batch = 5;


    int rank, numProcs;
    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);


    float ver;
    struct Zoltan_Struct *zz;
    MESH_DATA myTargets, mySources;

    mySources.numGlobalPoints = N * numProcs;
    mySources.numMyPoints = N;

    myTargets.numGlobalPoints = M * numProcs;
    myTargets.numMyPoints = M;

    if (Zoltan_Initialize(argc, argv, &ver) != ZOLTAN_OK) {
        printf("Zoltan failed to initialize. Exiting.\n");
        MPI_Finalize();
        exit(0);
    }


    zz = Zoltan_Create(MPI_COMM_WORLD);


    mySources.x = malloc(N*sizeof(double));
    mySources.y = malloc(N*sizeof(double));
    mySources.z = malloc(N*sizeof(double));
    double *qS = malloc(N*sizeof(double));
    double *wS = malloc(N*sizeof(double));

    myTargets.x = malloc(M*sizeof(double));
    myTargets.y = malloc(M*sizeof(double));
    myTargets.z = malloc(M*sizeof(double));
    double *qT = malloc(M*sizeof(double));

    double *potential = malloc(M*sizeof(double));

    for (int i = 0; i < N; ++i) {
        mySources.x[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        mySources.y[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        mySources.z[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        qS[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        wS[i] = 1.;
    }

    for (int i = 0; i < M; ++i) {
        myTargets.x[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        myTargets.y[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        myTargets.z[i] = ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        qT[i] = 1.;
    }

    treedriverWrapper(M, N, myTargets.x, myTargets.y, myTargets.z, qT,
                      mySources.x, mySources.y, mySources.z, qS, wS, potential,
                      pot_type, kappa, order, theta, max_per_leaf, max_per_batch);

    printf("Treedriver has finished.\n");

    free(mySources.x);
    free(mySources.y);
    free(mySources.z);
    free(qS);
    free(wS);
    free(myTargets.x);
    free(myTargets.y);
    free(myTargets.z);
    free(qT);
    free(potential);

    Zoltan_Destroy(&zz);
    MPI_Finalize();

    return 0;
}
