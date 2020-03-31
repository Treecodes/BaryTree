#include <mpi.h>


void START_TIMER(double *time)
{
    *time = MPI_Wtime();

    return;
}


void STOP_TIMER(double *time)
{
    *time = MPI_Wtime() - *time;

    return;
}
