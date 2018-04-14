#ifndef H_BATCH_H
#define H_BATCH_H

/* declaration of struct with tag batch */
struct batch
{
        int num;
        int *reorder;
        int **index;
        double **center;
        double *radius;
};

#endif /* H_BATCH_H */
