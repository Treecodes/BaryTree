#ifndef H_STRUCT_TREE_H
#define H_STRUCT_TREE_H

struct Tree
{
    int numnodes;
    int numleaves;
    
    int *ibeg;
    int *iend;
    int *numpar;
    
    int *cluster_ind;
    int *level;
    
    double *radius;

    double *x_mid;
    double *y_mid;
    double *z_mid;

    double *x_min;
    double *y_min;
    double *z_min;

    double *x_max;
    double *y_max;
    double *z_max;

    int *num_children;
    int *children;
};

#endif /* H_STRUCT_TREE_H */
