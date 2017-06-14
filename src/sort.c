/* sort functions for use by treecode routines */
#include <math.h>
#include "array.h"

void quicksortTargets(double *x, double *y, double *z, int *ind, int l, int r);
int partitionTargets(double *x, double *y, double *z, int *ind, int l, int r);

void sortTargets(double *x, double *y, double *z, int *ind, int numpars, int dflag)
{
    if (dflag == 0)
        quicksortTargets(x, y, z, ind, 0, numpars-1);
    else if (dflag == 1) 
        quicksortTargets(y, z, x, ind, 0, numpars-1);
    else
        quicksortTargets(z, x, y, ind, 0, numpars-1);

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


void interleaveGridTargets(double *x, double *y, double *z, int *ind, int numpars, int p)
{
    double *xx, *yy, *zz;
    int *indd;
    int modsize, moddiv;
    int numparsloc, maxparsdiv;
    int i;
    
    numparsloc = numpars/p;
    maxparsdiv = numpars/p * p;

    make_vector(xx, numpars);
    make_vector(yy, numpars);
    make_vector(zz, numpars);
    make_vector(indd, numpars);
    
    for (i = 0; i < maxparsdiv; i++) {
        modsize = i/p;
        moddiv = i%p;
        
        xx[numparsloc*moddiv + modsize] = x[i];
        yy[numparsloc*moddiv + modsize] = y[i];
        zz[numparsloc*moddiv + modsize] = z[i];
        indd[numparsloc*moddiv + modsize] = ind[i];
        
        //printf("%d : modsize %d, moddiv %d, newind %d\n", i, modsize, moddiv, numparsloc*moddiv + modsize);
    }
    
    for (i = 0; i < maxparsdiv; i++) {
        x[i] = xx[i];
        y[i] = yy[i];
        z[i] = zz[i];
        ind[i] = indd[i];
    }
    
    free_vector(xx);
    free_vector(yy);
    free_vector(zz);
    free_vector(indd);
    
    return;
}
