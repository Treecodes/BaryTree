/* sort functions for use by treecode routines */
#include <math.h>

void quicksortTargets(double *x, double *y, double *z, int *ind, int l, int r);
int partitionTargets(double *x, double *y, double *z, int *ind, int l, int r);

void sortTargets(double *x, double *y, double *z, int *ind, int numpars)
{
//    quicksortTargets(z, x, y, ind, 0, numpars-1);
//    quicksortTargets(y, x, z, ind, 0, numpars-1);
    quicksortTargets(x, y, z, ind, 0, numpars-1);

    return;
}


void quicksortTargets(double *x, double *y, double *z, int *ind, int l, int r)
{
    int i;

    if (l < r) {
        i = partitionTargets(x, y, z, ind, l, r);
        quicksortTargets(x, y, z, ind, l, i-1);
        quicksortTargets(x, y, z, ind, i+1, r);
    }

    return;
}


int partitionTargets(double *x, double *y, double *z, int *ind, int l, int r)
{
    double pivot, tx, ty, tz;
    int i, j, ti;

    pivot = x[l];
    i = l;
    j = r+1;

    while (1) {
        do ++i; while(x[i] <= pivot && i <= r);
        do --j; while(x[j] > pivot);
        if (i >= j) break;

          tx = x[i];   ty = y[i];   tz = z[i];     ti = ind[i];
        x[i] = x[j]; y[i] = y[j]; z[i] = z[j]; ind[i] = ind[j];
        x[j] = tx;   y[j] = ty;   z[j] = tz;   ind[j] = ti;
    }

      tx = x[l];   ty = y[l];   tz = z[l];     ti = ind[l];
    x[l] = x[j]; y[l] = y[j]; z[l] = z[j]; ind[l] = ind[j];
    x[j] = tx;   y[j] = ty;   z[j] = tz;   ind[j] = ti;

    return j;
}
