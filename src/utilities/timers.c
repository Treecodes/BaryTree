#include <mpi.h>


void start_timer(double *time)
{
    *time = MPI_Wtime();

    return;
}


void stop_timer(double *time)
{
    *time = MPI_Wtime() - *time;

    return;
}
