#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <zoltan.h>

#include "../src/treedriverWrapper.h"

int main(int argc, char **argv)
{

    int rank, numProcs;
    MPI_Init(&argc, &argv);

    int N = 100000;
    int M = 100000;

    int pot_type = 0; //Coulomb, Lagrange interpolation
    double kappa = 0.;

    int order = 5;
    double theta = 0.5;

    int max_per_leaf = 500;
    int max_per_batch = 5;

    int number_of_threads = 1;

    double *xS = malloc(N*sizeof(double));
    double *yS = malloc(N*sizeof(double));
    double *zS = malloc(N*sizeof(double));
    double *qS = malloc(N*sizeof(double));
    double *wS = malloc(N*sizeof(double));

    double *xT = malloc(M*sizeof(double));
    double *yT = malloc(M*sizeof(double));
    double *zT = malloc(M*sizeof(double));
    double *qT = malloc(M*sizeof(double));

    double *potential = malloc(M*sizeof(double));

    for (int i = 0; i < N; ++i) {
        xS[i] =  ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        yS[i] =  ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        zS[i] =  ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        qS[i] =  ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        wS[i] =  1.;
    }

    for (int i = 0; i < M; ++i) {
        xT[i] =  ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        yT[i] =  ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        zT[i] =  ((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        qT[i] =  1.;
    }

    treedriverWrapper(M, N, xT, yT, zT, qT, xS, yS, zS, qS, wS, potential,
                      pot_type, kappa, order, theta, max_per_leaf, max_per_batch);

    printf("Treedriver has finished.\n");

    free(xS);
    free(yS);
    free(zS);
    free(qS);
    free(wS);
    free(xT);
    free(yT);
    free(zT);
    free(qT);
    free(potential);

    MPI_Finalize();

    return 0;
}
