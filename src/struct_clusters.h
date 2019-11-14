#ifndef H_CLUSTERS_H
#define H_CLUSTERS_H

/* declaration of struct with tag particles */
struct clusters 
{
        int num;
        int num_weights;
        int num_charges;

        double *x;
        double *y;
        double *z;
        double *q;
        double *w;  // quadrature weights.  Set equal to 1 if interacting with particles, not performing convolution integral.
};

#endif /* H_CLUSTERS_H */
