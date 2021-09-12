#ifndef H_CLUSTERS_H
#define H_CLUSTERS_H

/* declaration of struct with tag particles */
struct Clusters
{
        int num;
        int num_charges;

        double *x;
        double *y;
        double *z;
        double *q;
};

#endif /* H_CLUSTERS_H */
