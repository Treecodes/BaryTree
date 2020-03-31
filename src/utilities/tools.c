/* tool functions for use by treecode routines */
#include <math.h>
#include "array.h"


double minval(double *x, int numels)
{
    double min = x[0];

    for (int i = 1; i < numels; ++i) {
        if (min > x[i])
            min = x[i];
    }

    return min;
}



double maxval(double *x, int numels)
{
    double max = x[0];

    for (int i = 1; i < numels; ++i) {
        if (max < x[i])
            max = x[i];
    }

    return max;
}



double sum(double *x, int numels)
{
    double sum = 0.0;

    for (int i = 0; i < numels; ++i)
        sum += x[i];

    return sum;
}



int sum_int(int *x, int numels)
{
    int sum = 0.0;

    for (int i = 0; i < numels; ++i)
        sum += x[i];

    return sum;
}



double max3(double a, double b, double c)
{
    double max = a;

    if (max < b) max = b;
    if (max < c) max = c;

    return max;
}



double min3(double a, double b, double c)
{
    double min = a;

    if (min > b) min = b;
    if (min > c) min = c;

    return min;
}



int maxval_int(int *x, int numels)
{
    int max = x[0];

    for (int i = 1; i < numels; i++) {
        if (max < x[i])
            max = x[i];
    }

    return max;
}
