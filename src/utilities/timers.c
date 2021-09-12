#include <omp.h>


void START_TIMER(double *time)
{
    *time = omp_get_wtime();

    return;
}


void STOP_TIMER(double *time)
{
    *time = omp_get_wtime() - *time;

    return;
}
