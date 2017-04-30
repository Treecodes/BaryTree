/* sort functions for use by treecode routines */
#include <math.h>

void quicksortTargets(double *x, double *y, double *z, int l, int r);
int partitionTargets(double *x, double *y, double *z, int l, int r);

void sortTargets(double *x, double *y, double *z, int numpars)
{
//    quicksortTargets(z, x, y, 0, numpars-1);
//    quicksortTargets(y, x, z, 0, numpars-1);
    quicksortTargets(x, y, z, 0, numpars-1);

    return;
}


void quicksortTargets(double *x, double *y, double *z, int l, int r)
{
    int i;

    if (l < r)
    {
        i = partitionTargets(x, y, z, l, r);
        quicksortTargets(x, y, z, l, i-1);
        quicksortTargets(x, y, z, i+1, r);
    }

    return;
}


int partitionTargets(double *x, double *y, double *z, int l, int r)
{
    double pivot, tx, ty, tz;
    int i, j;

    pivot = x[l];
    i = l;
    j = r+1;

    while (1)
    {
        do ++i; while(x[i] <= pivot && i <= r);
        do --j; while(x[j] > pivot);
        if (i >= j) break;

          tx = x[i];   ty = y[i];   tz = z[i];
        x[i] = x[j]; y[i] = y[j]; z[i] = z[j];
        x[j] = tx;   y[j] = ty;   z[j] = tz;
    }

      tx = x[l];   ty = y[l];   tz = z[l];
    x[l] = x[j]; y[l] = y[j]; z[l] = z[j];
    x[j] = tx;   y[j] = ty;   z[j] = tz;

    return j;
}
