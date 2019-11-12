/* file minunit_example.c */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <zoltan.h>
#include <time.h>
#include <float.h>

#include "minunit.h"
#include "../src/treedriver.h"
#include "../src/particles.h"


int tests_run = 0;

int foo = 7;
int bar = 5;

static char * test_foo() {
    mu_assert("error, foo != 7", foo == 7);
    return 0;
}

static char * test_bar() {
    mu_assert("error, bar != 5", bar == 5);
    return 0;
}

static char * test_direct_sum_on_10_particles() {

    int N=10;

    struct particles *sources = NULL;
    struct particles *targets = NULL;
    int *particleOrder = NULL;
    double *potential = NULL;
    double potential_engy = 0;

    sources = malloc(sizeof(struct particles));
    targets = malloc(sizeof(struct particles));
    potential = malloc(sizeof(double) * N);
    particleOrder = malloc(sizeof(int) * N);

    targets->num = N;
    targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));
    targets->order = malloc(targets->num*sizeof(int));

    sources->num = N;
    sources->x = malloc(sources->num*sizeof(double));
    sources->y = malloc(sources->num*sizeof(double));
    sources->z = malloc(sources->num*sizeof(double));
    sources->q = malloc(sources->num*sizeof(double));
    sources->w = malloc(sources->num*sizeof(double));
    sources->order = malloc(sources->num*sizeof(int));


    for (int i=0; i<targets->num; i++){

        targets->x[i]=1.0*i;
        targets->y[i]=1.0*i;
        targets->z[i]=1.0*i;
        targets->q[i]=1.0*i;
        targets->order[i] = i;

        sources->x[i]=1.0*i;
        sources->y[i]=1.0*i;
        sources->z[i]=1.0*i;
        sources->q[i]=1.0*i;
        sources->w[i]=1.0;
        sources->order[i] = i;

        potential[i]=0.0;
        }

    int max_per_leaf=100;
    int max_per_batch=100;
    double time_tree[9];

    char *kernelName, *singularityHandling, *approximationName;

    kernelName="coulomb";
    singularityHandling="skipping";
    approximationName="lagrange";
    int tree_type=1;
    double kappa=0.5;

    treedriver(sources, targets, 0, 0.0, max_per_leaf, max_per_batch,
                   kernelName, kappa, singularityHandling, approximationName, tree_type, potential,
                   &potential_engy, time_tree);

    for (int i=0; i<targets->num; i++){
        double trueValue=0.0;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += j/(r);
        }
        for (int j=i+1; j<targets->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += j/(r);
        }

        mu_assert("Direct sum potential not correct for coulomb kernel with skipping", abs(potential[i] - trueValue) < 1e-10);
    }


    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    singularityHandling="subtraction";
    treedriver(sources, targets, 0, 0.0, max_per_leaf, max_per_batch,
                   kernelName, kappa, singularityHandling, approximationName, tree_type, potential,
                   &potential_engy, time_tree);

    for (int i=0; i<targets->num; i++){
        double trueValue=0.0;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa*kappa) )/(r);
        }
        for (int j=i+1; j<targets->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa*kappa) )/(r);
        }

        mu_assert("Direct sum potential not correct for coulomb kernel with subtraction", abs(potential[i] - trueValue) < 1e-10);
    }

    kernelName="yukawa";
    singularityHandling="skipping";

    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    treedriver(sources, targets, 0, 0.0, max_per_leaf, max_per_batch,
                   kernelName, kappa, singularityHandling, approximationName, tree_type, potential,
                   &potential_engy, time_tree);

    for (int i=0; i<targets->num; i++){
        double trueValue=0.0;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += j*exp(-kappa*r)/(r);
        }
        for (int j=i+1; j<targets->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += j*exp(-kappa*r)/(r);
        }

        mu_assert("Direct sum potential not correct for yukawa kernel with skipping", abs(potential[i] - trueValue) < 1e-10);
    }


    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    singularityHandling="subtraction";
    treedriver(sources, targets, 0, 0.0, max_per_leaf, max_per_batch,
                   kernelName, kappa, singularityHandling, approximationName, tree_type, potential,
                   &potential_engy, time_tree);

    for (int i=0; i<targets->num; i++){
        double trueValue=0.0;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }
        for (int j=i+1; j<targets->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }

        mu_assert("Direct sum potential not correct for yukawa kernel with subtraction", abs(potential[i] - trueValue) < 1e-10);
    }



    free(sources->x);
    free(sources->y);
    free(sources->z);
    free(sources->q);
    free(sources->w);
    free(sources->order);
    free(sources);

    free(targets->x);
    free(targets->y);
    free(targets->z);
    free(targets->q);
    free(targets->order);
    free(targets);

    free(potential);
    return 0;
}

static char * all_tests() {
    mu_run_test(test_foo);
    mu_run_test(test_bar);
    mu_run_test(test_direct_sum_on_10_particles);
return 0;
}

int main(int argc, char **argv) {
    int rc, rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);


    char *result = all_tests();
    if (result != 0) {
        printf("%s\n", result);
    }
    else {
        printf("ALL TESTS PASSED\n");
    }
        printf("Tests run: %d\n", tests_run);


    MPI_Finalize();
    return result != 0;
}
