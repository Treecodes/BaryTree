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
    mu_assert("TEST FAILED: error, foo != 7", foo == 7);
    return 0;
}

static char * test_bar() {
    mu_assert("TEST FAILED: error, bar != 5", bar == 5);
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

        mu_assert("TEST FAILED: Direct sum potential not correct for coulomb kernel with skipping", abs(potential[i] - trueValue) < 1e-10);
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

        mu_assert("TEST FAILED: Direct sum potential not correct for coulomb kernel with subtraction", abs(potential[i] - trueValue) < 1e-10);
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

        mu_assert("TEST FAILED: Direct sum potential not correct for yukawa kernel with skipping", abs(potential[i] - trueValue) < 1e-10);
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

        mu_assert("TEST FAILED: Direct sum potential not correct for yukawa kernel with subtraction", abs(potential[i] - trueValue) < 1e-10);
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

static char * test_treecode_on_100_particles() {

    int N=100;

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

    int max_per_leaf=3;
    int max_per_batch=3;
    double time_tree[9];

    char *kernelName, *singularityHandling, *approximationName;

    kernelName="coulomb";
    singularityHandling="skipping";
    int tree_type=1;
    double kappa=0.5;

    int order=2;
    double theta=0.7;

    //////////// Lagrange /////////
    approximationName="lagrange";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
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
//        printf("potential at %i:, %f\n", i, potential[i]);
//        printf("trueValue at %i:, %f\n\n", i, trueValue);
        mu_assert("TEST FAILED: TEST FAILED: Treecode potential not correct for: lagrange-coulomb-skipping", fabs(potential[i] - trueValue)/fabs(trueValue) < 3e-3);
    }


    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    singularityHandling="subtraction";
    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
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

//        printf("potential at %i:, %f\n", i, potential[i]);
//        printf("trueValue at %i:, %f\n", i, trueValue);
//        printf("Absolute error: %f\n", fabs( (potential[i] - trueValue) ) );
//        printf("Relative error at %i:: %f\n", i, fabs( (potential[i] - trueValue)/trueValue ) );
//        mu_assert("TEST FAILED: TEST FAILED: Treecode potential not correct for coulomb kernel with subtraction", (abs(potential[i] - trueValue)/abs(trueValue)) < 1e-100);
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-coulomb-subtraction", fabs(potential[i] - trueValue)/fabs(trueValue) < 2e-2);
    }

    kernelName="yukawa";
    singularityHandling="skipping";

    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
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

//        mu_assert("TEST FAILED: Direct sum potential not correct for yukawa kernel with skipping", abs(potential[i] - trueValue) < 1e-10);
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-yukawa-skipping", fabs(potential[i] - trueValue)/fabs(trueValue) < 8e-3);
    }


    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    singularityHandling="subtraction";
    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
                   kernelName, kappa, singularityHandling, approximationName, tree_type, potential,
                   &potential_engy, time_tree);

    for (int i=0; i<targets->num; i++){
        double trueValue=0.0;
        for (int j=0; j<i; j++){
            double r = fabs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }
        for (int j=i+1; j<sources->num; j++){
            double r = fabs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }
//        printf("potential at %i:, %1.2e\n", i, potential[i]);
//        printf("trueValue at %i:, %1.2e\n", i, trueValue);
//        printf("Absolute error at %i: %f\n\n", i, fabs( (potential[i] - trueValue) ) );
//        printf("Relative error at %i: %f\n\n", i, fabs( (potential[i] - trueValue)/trueValue ) );

        // measure absolute error for this example, since true values are very close to zero.
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-yukawa-subtraction", fabs(potential[i] - trueValue) < 2e-2);

    }


    //////// Hermite ////////

    approximationName="hermite";
    kernelName="coulomb";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
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
//        printf("potential at %i:, %f\n", i, potential[i]);
//        printf("trueValue at %i:, %f\n\n", i, trueValue);
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-coulomb-skipping", fabs(potential[i] - trueValue)/fabs(trueValue) < 3e-4);
    }



    approximationName="hermite";
    kernelName="coulomb";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
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

//        printf("potential at %i:, %f\n", i, potential[i]);
//        printf("trueValue at %i:, %f\n", i, trueValue);
//        printf("Absolute error: %f\n", fabs( (potential[i] - trueValue) ) );
//        printf("Relative error at %i:: %f\n", i, fabs( (potential[i] - trueValue)/trueValue ) );
//        mu_assert("TEST FAILED: Treecode potential not correct for coulomb kernel with subtraction", (abs(potential[i] - trueValue)/abs(trueValue)) < 1e-100);
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-coulomb-subtraction", fabs(potential[i] - trueValue)/fabs(trueValue) < 2e-2);
    }

    approximationName="hermite";
    kernelName="yukawa";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
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

//        mu_assert("TEST FAILED: Direct sum potential not correct for yukawa kernel with skipping", abs(potential[i] - trueValue) < 1e-10);
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-yukawa-skipping", fabs(potential[i] - trueValue)/fabs(trueValue) < 5e-4);
    }


    approximationName="hermite";
    kernelName="yukawa";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
                   kernelName, kappa, singularityHandling, approximationName, tree_type, potential,
                   &potential_engy, time_tree);

    for (int i=0; i<targets->num; i++){
        double trueValue=0.0;
        for (int j=0; j<i; j++){
            double r = fabs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }
        for (int j=i+1; j<sources->num; j++){
            double r = fabs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }
//        printf("potential at %i:, %1.2e\n", i, potential[i]);
//        printf("trueValue at %i:, %1.2e\n", i, trueValue);
//        printf("Absolute error at %i: %f\n\n", i, fabs( (potential[i] - trueValue) ) );
//        printf("Relative error at %i: %f\n\n", i, fabs( (potential[i] - trueValue)/trueValue ) );

        // measure absolute error for this example, since true values are very close to zero.
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-yukawa-subtraction", fabs(potential[i] - trueValue) < 2e-3);

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
    mu_run_test(test_treecode_on_100_particles);
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
