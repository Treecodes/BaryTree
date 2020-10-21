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

};

#endif /* H_STRUCT_TREE_H */
