#ifndef H_BATCH_H
#define H_BATCH_H

/* declaration of struct with tag batch */
struct batch
{
        //int num;
        int numnodes;

        int *reorder;

        // int **index;
        int *ibeg;
        int *iend;
        int *numpar;

        int *numApprox;
        int *numDirect;

        // double **center;
        double *x_mid;
        double *y_mid;
        double *z_mid;

        double *radius;

        int *numClusterCharges;
        int *numClusterWeights;
};

#endif /* H_BATCH_H */
