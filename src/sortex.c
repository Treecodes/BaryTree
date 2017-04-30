#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "sort.h"

int main(int argc, char **argv)
{
    int rank, p;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    double *x, *y, *z;
    int i, numpars;
    int numparsloc, maxparsloc, globparsloc;

    numpars = atoi(argv[1]);
    
    numparsloc = (int)floor((double)numpars/(double)p);
    maxparsloc = numparsloc + (numpars - (int)floor((double)numpars/(double)p) * p);
    
    if (rank == 0) numparsloc = maxparsloc;
    globparsloc = maxparsloc + numparsloc * (rank-1);
    
    
    //Construct unsorted list
    x = malloc(maxparsloc * sizeof(double));
    y = malloc(maxparsloc * sizeof(double));
    z = malloc(maxparsloc * sizeof(double));

    for (i = 0; i < numparsloc; i++)
    {
        x[i] = (double)(rand() % 5);
        y[i] = (double)(rand() % 5);
        z[i] = (double)(rand() % 5);
    }

    
    
    
    
    sortTargets(x, y, z, numparsloc);

    
    
    
    
    //Print list
    for (i = 0; i < numparsloc; i++)
    {
        printf("proc %d: %f,  %f,  %f\n", rank, x[i], y[i], z[i]);
    }

    MPI_Finalize();
    
    return 0;
}
