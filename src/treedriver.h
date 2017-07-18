#ifndef H_TREEDRIVER_H
#define H_TREEDRIVER_H

/* declaration of primary treecode driver */

void treecode(double *xS, double *yS, double *zS, double *qS, 
              double *xT, double *yT, double *zT,
              int numparsS, int numparsT, double *tEn, double *tpeng, 
              int order, double theta, int shrink, int maxparnode,
              double timetree[4], int treelevel, int iflag,
              int pot_type, double kappa, int tree_type);


void treecode_grid(double *xS, double *yS, double *zS, double *qS, 
                   double *xyzminmax, int *xyzdim, int numparsS,
                   double *tEn, double *tpeng, int order, double theta, 
                   int maxparnode, double *timetree, int pot_type, double kappa);


void treecode_grid_bdry(double *xS, double *yS, double *zS, double *qS,
                   double *xyzminmax, int *xyzdim, double zyx, int dir, int numparsS,
                   double *tEn, double *tpeng, int order, double theta,
                   int maxparnode, double *timetree, int pot_type, double kappa);

#endif /* H_TREEFUNCTIONS_H */
