#ifndef H_TREEDRIVER_H
#define H_TREEDRIVER_H

/* declaration of primary treecode driver */

void treecode(double *xS, double *yS, double *zS, double *qS, 
              double *xT, double *yT, double *zT,
              int numparsS, int numparsT, double *tEn, double *tpeng, 
              int order, double theta, int shrinkS, int shrinkT, 
              int maxparnodeS, int maxparnodeT, double *timetree,
              int treelevelS, int treelevelT, int iflagS, int iflagT,
              int pot_type, double kappa);

#endif /* H_TREEFUNCTIONS_H */
