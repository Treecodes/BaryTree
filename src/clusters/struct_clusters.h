#ifndef H_CLUSTERS_H
#define H_CLUSTERS_H

/* declaration of struct with tag particles */
struct Clusters
{
        int num;
        int num_weights;
        int num_charges;

        double *x;
        double *y;
        double *z;
        double *q;
        // quadrature weights.  Set = 1 if interacting particles, not performing convolution integral.
        double *w;
};

#endif /* H_CLUSTERS_H */
