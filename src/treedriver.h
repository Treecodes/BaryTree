#ifndef H_TREEDRIVER_H
#define H_TREEDRIVER_H

/* declaration of primary treecode driver */

void treedriver(double *xS, double *yS, double *zS, double *qS,
                double *xT, double *yT, double *zT,
                int numparsS, int numparsT,
                int order, double theta, int maxparnode, int batch_size,
                int pot_type, double kappa, int tree_type,
                double *tEn, double *tpeng, double *timetree);

#endif /* H_TREEFUNCTIONS_H */
