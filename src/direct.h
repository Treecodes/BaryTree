#ifndef H_DIRECT_H
#define H_DIRECT_H

void direct_eng_f(float *xS, float *yS, float *zS, float *qS, float *wS,
                float *xT, float *yT, float *zT, float *qT,
                int numparsS, int numparsT, float *denergy, float *dpeng,
                int pot_type, float kappa);


void direct_eng_d(double *xS, double *yS, double *zS, double *qS, double *wS,
                double *xT, double *yT, double *zT, double *qT,
                int numparsS, int numparsT, double *denergy, double *dpeng,
                int pot_type, double kappa);

#endif /* H_DIRECT_H */
