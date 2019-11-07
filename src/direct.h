#ifndef H_DIRECT_H
#define H_DIRECT_H

void direct_eng(double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                double kappa, char *kernelName);

#endif /* H_DIRECT_H */
