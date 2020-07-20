#ifndef H_STRUCT_TREE_H
#define H_STRUCT_TREE_H

struct Tree
{
    int numnodes;
    int numleaves;
    
    int min_leaf_size;
    int max_leaf_size;
    int max_depth;
    
    int *ibeg;
    int *iend;
    int *numpar;
    
    int *cluster_ind;
    int *used;
    int *used_children;
    int *used_parent;
    int *used_leaf;
    
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
    int *parent;

    int **levels_list;
    int *levels_list_num;

    int *leaves_list;
    int leaves_list_num;

    int *x_dim;
    int *y_dim;
    int *z_dim;

    int *x_low_ind;
    int *y_low_ind;
    int *z_low_ind;

    int *x_high_ind;
    int *y_high_ind;
    int *z_high_ind;
};

#endif /* H_STRUCT_TREE_H */
