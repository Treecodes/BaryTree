/* tool functions for use by treecode routines */
#ifndef H_TOOLS_H                                                                     
#define H_TOOLS_H 

double minval(double *x, int numels);
double maxval(double *x, int numels);
int maxval_int(int *x, int numels);

double sum(double *x, int numels);
int sum_int(int *x, int numels);

double max3(double a, double b, double c);
double min3(double a, double b, double c);

#endif /* H_TOOLS_H */
