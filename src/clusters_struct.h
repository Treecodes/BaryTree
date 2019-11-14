#ifndef H_TNODE_H
#define H_TNODE_H

/* declaration of struct with tag tnode */
struct clusters
{
    int numnodes;
    
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



    int *numApprox;
    int *numDirect;
    int *reorder;
};


#endif /* H_TNODE_H */
