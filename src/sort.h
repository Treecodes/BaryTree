/* sort functions for use on targets */
#ifndef H_SORT_H                                                                     
#define H_SORT_H 

void sortTargets(double *xT, double *yT, double *zT, double *qT, int *iT, int numparsT, int dflag);
void interleaveGridTargets(double *x, double *y, double *z, double *q, int *ind, int numpars, int p);

#endif /* H_SORT_H */
