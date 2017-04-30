#ifndef H_TREEDRIVER_H
#define H_TREEDRIVER_H

/* declaration of primary treecode driver */

void treecode(double *xS, double *yS, double *zS, double *qS, 
              double *xT, double *yT, double *zT,
              int numparsS, int numparsT, double *tEn, double *tpeng, 
              int order, double theta, int shrink, int maxparnode,
              double *timetree, int treelevel, int iflag,
              int pot_type, double kappa, int tree_type);

#endif /* H_TREEFUNCTIONS_H */
