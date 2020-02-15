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
#include "../src/directdriver.h"
#include "../src/struct_particles.h"
#include "../src/struct_kernel.h"
#include "../src/kernel.h"
#include "../src/array.h"


int tests_run = 0;

static char * test_direct_sum_on_10_particles() {

    int N=10;

    int verbosity=0;

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

//    kernelName="coulomb";

    struct kernel *kernel = NULL;
    kernel = malloc(sizeof (struct kernel));
    double * kernelParameters = NULL;
    int numberOfKernelParameters=1;
    make_vector(kernelParameters,numberOfKernelParameters);
    kernelParameters[0]=0.5;
    AllocateKernelStruct(kernel, numberOfKernelParameters, "coulomb");

    int tree_type=1;
    double kappa=0.5;
    SetKernelParameters(kernel, kernelParameters);
    singularityHandling="skipping";

    directdriver(sources, targets, kernel, singularityHandling, approximationName,
                 potential, time_tree);

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

        mu_assert("TEST FAILED: Direct sum potential not correct for coulomb kernel with skipping", fabs(potential[i] - trueValue) < 1e-10);
    }


    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    singularityHandling="subtraction";
    kernelParameters[0]=5.5;
    kappa=5.5;
    SetKernelParameters(kernel, kernelParameters);
    directdriver(sources, targets, kernel, singularityHandling, approximationName,
                 potential, time_tree);

    for (int i=0; i<targets->num; i++){
        double trueValue=2.0 * M_PI * kappa * kappa * i;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa/kappa) )/(r);
        }
        for (int j=i+1; j<targets->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa/kappa) )/(r);
        }

        mu_assert("TEST FAILED: Direct sum potential not correct for coulomb kernel with subtraction", fabs(potential[i] - trueValue) < 1e-10);
    }

    singularityHandling="skipping";
    struct kernel *kernel2 = NULL;
    kernel2 = malloc(sizeof (struct kernel));
    numberOfKernelParameters=1;
    kappa=0.5;
    AllocateKernelStruct(kernel2, numberOfKernelParameters, "yukawa");
    SetKernelParameters(kernel2, &kappa);

    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }

    directdriver(sources, targets, kernel2, singularityHandling, approximationName,
                 potential, time_tree);

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

        mu_assert("TEST FAILED: Direct sum potential not correct for yukawa kernel with skipping", fabs(potential[i] - trueValue) < 1e-10);
    }


    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    singularityHandling="subtraction";
    directdriver(sources, targets, kernel2, singularityHandling, approximationName,
                 potential, time_tree);

    for (int i=0; i<targets->num; i++){
        double trueValue=4.0 * M_PI / kappa / kappa * i;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }
        for (int j=i+1; j<targets->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }

        mu_assert("TEST FAILED: Direct sum potential not correct for yukawa kernel with subtraction", fabs(potential[i] - trueValue) < 1e-10);
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

    FreeKernelStruct(kernel);
    FreeKernelStruct(kernel2);
    return 0;
}

static char * test_treecode_on_100_particles() {


    struct kernel *CoulombKernel = NULL;
    struct kernel *YukawaKernel = NULL;
    CoulombKernel = malloc(sizeof (struct kernel));
    YukawaKernel = malloc(sizeof (struct kernel));
    double * kernelParameters = NULL;
    int numberOfKernelParameters=1;
    make_vector(kernelParameters,numberOfKernelParameters);
    AllocateKernelStruct(CoulombKernel, numberOfKernelParameters, "coulomb");
    AllocateKernelStruct(YukawaKernel, numberOfKernelParameters, "yukawa");


    int verbosity=0;
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
    SetKernelParameters(CoulombKernel, &kappa);
    SetKernelParameters(YukawaKernel, &kappa);

    int order=2;
    double theta=0.7;

    /***********************************************/
    /******************* Test 1 ********************/
    /***********************************************/
    /********** lagrange-coulomb-skipping **********/
    /***********************************************/
    approximationName="lagrange";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

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

        if (verbosity>0) printf("trueValue = %1.8e\n", trueValue);
        if (verbosity>0) printf("computedValue = %1.8e\n", potential[i]);
        if (verbosity>0) printf("relative error = %1.8e\n\n", fabs(potential[i] - trueValue)/fabs(trueValue));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-coulomb-skipping", fabs(potential[i] - trueValue)/fabs(trueValue) < 3e-3);
    }


    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }

    /***********************************************/
    /******************* Test 2 ********************/
    /***********************************************/
    /********* lagrange-coulomb-subtraction ********/
    /***********************************************/
    singularityHandling="subtraction";

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        double trueValue=2.0 * M_PI * kappa * kappa * i;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa*kappa) )/(r);
        }
        for (int j=i+1; j<targets->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa*kappa) )/(r);
        }
        if (verbosity>0) printf("trueValue = %1.8e\n", trueValue);
        if (verbosity>0) printf("computedValue = %1.8e\n", potential[i]);
        if (verbosity>0) printf("relative error = %1.8e\n\n", fabs(potential[i] - trueValue)/fabs(trueValue));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-coulomb-subtraction", fabs(potential[i] - trueValue)/fabs(trueValue) < 2e-2);
    }

    /***********************************************/
    /******************* Test 3 ********************/
    /***********************************************/
    /*********** lagrange-yukawa-skipping **********/
    /***********************************************/
    kernelName="yukawa";
    singularityHandling="skipping";

    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

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
        if (verbosity>0) printf("trueValue = %1.8e\n", trueValue);
        if (verbosity>0) printf("computedValue = %1.8e\n", potential[i]);
        if (verbosity>0) printf("relative error = %1.8e\n\n", fabs(potential[i] - trueValue)/fabs(trueValue));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-yukawa-skipping", fabs(potential[i] - trueValue)/fabs(trueValue) < 8e-3);
    }

    /***********************************************/
    /******************* Test 4 ********************/
    /***********************************************/
    /********* lagrange-yukawa-subtraction *********/
    /***********************************************/
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    singularityHandling="subtraction";

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        double trueValue=4.0 * M_PI / kappa / kappa * i;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }
        for (int j=i+1; j<sources->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }
        // measure absolute error for this example, since true values are very close to zero.
        if (verbosity>0) printf("trueValue = %1.8e\n", trueValue);
        if (verbosity>0) printf("computedValue = %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error = %1.8e\n\n", fabs(potential[i] - trueValue));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-yukawa-subtraction", fabs(potential[i] - trueValue) < 2e-2);

    }

    /***********************************************/
    /******************* Test 5 ********************/
    /***********************************************/
    /********* hermite-coulomb-skipping ************/
    /***********************************************/
    approximationName="hermite";
    kernelName="coulomb";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

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
        if (verbosity>0) printf("trueValue = %1.8e\n", trueValue);
        if (verbosity>0) printf("computedValue = %1.8e\n", potential[i]);
        if (verbosity>0) printf("relative error = %1.8e\n\n", fabs(potential[i] - trueValue)/fabs(trueValue));
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-coulomb-skipping", fabs(potential[i] - trueValue)/fabs(trueValue) < 3e-4);
    }


    /***********************************************/
    /******************* Test 6 ********************/
    /***********************************************/
    /******** hermite-coulomb-subtraction **********/
    /***********************************************/
    approximationName="hermite";
    kernelName="coulomb";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        double trueValue=2.0 * M_PI * kappa * kappa * i;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa*kappa) )/(r);
        }
        for (int j=i+1; j<targets->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa*kappa) )/(r);
        }
        if (verbosity>0) printf("trueValue = %1.8e\n", trueValue);
        if (verbosity>0) printf("computedValue = %1.8e\n", potential[i]);
        if (verbosity>0) printf("relative error = %1.8e\n\n", fabs(potential[i] - trueValue)/fabs(trueValue));
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-coulomb-subtraction", fabs(potential[i] - trueValue)/fabs(trueValue) < 2e-2);
    }

    /***********************************************/
    /******************* Test 7 ********************/
    /***********************************************/
    /********** hermite-yukawa-skipping ************/
    /***********************************************/
    approximationName="hermite";
    kernelName="yukawa";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

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
        if (verbosity>0) printf("trueValue = %1.8e\n", trueValue);
        if (verbosity>0) printf("computedValue = %1.8e\n", potential[i]);
        if (verbosity>0) printf("relative error = %1.8e\n\n", fabs(potential[i] - trueValue)/fabs(trueValue));
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-yukawa-skipping", fabs(potential[i] - trueValue)/fabs(trueValue) < 5e-4);
    }

    /***********************************************/
    /******************* Test 8 ********************/
    /***********************************************/
    /********* hermite-yukawa-subtraction **********/
    /***********************************************/
    approximationName="hermite";
    kernelName="yukawa";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
    }
    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        double trueValue=4.0 * M_PI / kappa / kappa * i;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }
        for (int j=i+1; j<sources->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i)*exp(-kappa*r)/(r);
        }
        // measure absolute error for this example, since true values are very close to zero.
        if (verbosity>0) printf("trueValue = %1.8e\n", trueValue);
        if (verbosity>0) printf("computedValue = %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error = %1.8e\n\n", fabs(potential[i] - trueValue));
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

    FreeKernelStruct(CoulombKernel);
    FreeKernelStruct(YukawaKernel);
    return 0;
}


static char * test_treecode_on_1_target_10000_sources() {

    int N=10000;
    int verbosity=0;

    struct kernel *CoulombKernel = NULL;
    struct kernel *YukawaKernel = NULL;
    CoulombKernel = malloc(sizeof (struct kernel));
    YukawaKernel = malloc(sizeof (struct kernel));
    double * kernelParameters = NULL;
    int numberOfKernelParameters=1;
    make_vector(kernelParameters,numberOfKernelParameters);
    AllocateKernelStruct(CoulombKernel, numberOfKernelParameters, "coulomb");
    AllocateKernelStruct(YukawaKernel, numberOfKernelParameters, "yukawa");


    struct particles *sources = NULL;
    struct particles *targets = NULL;
    int *particleOrder = NULL;
    double *potential = NULL, *potential_direct = NULL;
    double potential_engy = 0;
    double potential_engy_direct = 0;

    sources = malloc(sizeof(struct particles));
    targets = malloc(sizeof(struct particles));
    potential = malloc(sizeof(double) * N);
    potential_direct = malloc(sizeof(double) * N);
    particleOrder = malloc(sizeof(int) * N);

    targets->num = 1; //single target
    targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));
    targets->order = malloc(targets->num*sizeof(int));


    sources->num = N; // 10,000 sources
    sources->x = malloc(sources->num*sizeof(double));
    sources->y = malloc(sources->num*sizeof(double));
    sources->z = malloc(sources->num*sizeof(double));
    sources->q = malloc(sources->num*sizeof(double));
    sources->w = malloc(sources->num*sizeof(double));
    sources->order = malloc(sources->num*sizeof(int));


    for (int i=0; i<targets->num; i++){
        // single target at origin
        targets->x[i]=0.0;
        targets->y[i]=0.0;
        targets->z[i]=0.0;
        targets->q[i]=1.0;
        targets->order[i] = i;
        }

    srand(1);
    for (int i=0; i<sources->num; i++){
        // 10,000 randomly distributed sources in the [-1,1] box
        sources->x[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->y[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->z[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->q[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->w[i]=((double)rand()/(double)(RAND_MAX));
        sources->order[i] = i;
        }


    int max_per_leaf=100;
    int max_per_batch=100;
    double time_tree[9];

    char *kernelName, *singularityHandling, *approximationName;


    int tree_type=1; // particle-cluster
    double kappa=0.5;
    kernelParameters[0]=kappa;
    SetKernelParameters(CoulombKernel, &kappa);
    SetKernelParameters(YukawaKernel, &kappa);

    int order=4;
    double theta=0.8;

    /***********************************************/
    /******************* Test 1 ********************/
    /***********************************************/
    /********** lagrange-coulomb-skipping **********/
    /***********************************************/
    approximationName="lagrange";
    kernelName="coulomb";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, CoulombKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nlagrange-coulomb-skipping\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-coulomb-skipping", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-4);
    }



    /***********************************************/
    /******************* Test 2 ********************/
    /***********************************************/
    /********* lagrange-coulomb-subtraction ********/
    /***********************************************/
    approximationName="lagrange";
    kernelName="coulomb";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, CoulombKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nlagrange-coulomb-subtraction\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-coulomb-subtraction", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-5);
    }

    /***********************************************/
    /******************* Test 3 ********************/
    /***********************************************/
    /*********** lagrange-yukawa-skipping **********/
    /***********************************************/
    approximationName="lagrange";
    kernelName="yukawa";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, YukawaKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nlagrange-yukawa-skipping\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-yukawa-skipping", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-4);
    }

    /***********************************************/
    /******************* Test 4 ********************/
    /***********************************************/
    /********* lagrange-yukawa-subtraction *********/
    /***********************************************/
    approximationName="lagrange";
    kernelName="yukawa";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, YukawaKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nlagrange-yukawa-subtraction\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-yukawa-subtraction", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 6e-6);
    }

    /***********************************************/
    /******************* Test 5 ********************/
    /***********************************************/
    /********* hermite-coulomb-skipping ************/
    /***********************************************/
    approximationName="hermite";
    kernelName="coulomb";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, CoulombKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);
    if (verbosity>0) printf("\nhermite-coulomb-skipping finished direct reference.\n");

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);
    if (verbosity>0) printf("\nhermite-coulomb-skipping finished treecode run.\n");

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nhermite-coulomb-skipping\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-coulomb-skipping", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 3e-8);
    }


    /***********************************************/
    /******************* Test 6 ********************/
    /***********************************************/
    /******** hermite-coulomb-subtraction **********/
    /***********************************************/
    approximationName="hermite";
    kernelName="coulomb";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, CoulombKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nhermite-coulomb-subtraction\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-coulomb-subtraction", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-7);
    }

    /***********************************************/
    /******************* Test 7 ********************/
    /***********************************************/
    /********** hermite-yukawa-skipping ************/
    /***********************************************/
    approximationName="hermite";
    kernelName="yukawa";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, YukawaKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nhermite-yukawa-skipping\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]));
        mu_assert("TEST FAILED: TEST FAILED: Treecode potential not correct for: hermite-yukawa-skipping", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 4e-8);
    }

    /***********************************************/
    /******************* Test 8 ********************/
    /***********************************************/
    /********* hermite-yukawa-subtraction **********/
    /***********************************************/
    approximationName="hermite";
    kernelName="yukawa";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, YukawaKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order, theta, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential, time_tree, 1.0, verbosity);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nhermite-yukawa-subtraction\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]));
        mu_assert("TEST FAILED: TEST FAILED: Treecode potential not correct for: hermite-yukawa-subtraction", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 3e-8);
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
    free(potential_direct);
    return 0;
}


static char * test_treecode_parameters_on_1_target_10000_sources() {

    struct kernel *CoulombKernel = NULL;
    struct kernel *YukawaKernel = NULL;
    CoulombKernel = malloc(sizeof (struct kernel));
    YukawaKernel = malloc(sizeof (struct kernel));
    double * kernelParameters = NULL;
    int numberOfKernelParameters=1;
    make_vector(kernelParameters,numberOfKernelParameters);
    AllocateKernelStruct(CoulombKernel, numberOfKernelParameters, "coulomb");
    AllocateKernelStruct(YukawaKernel, numberOfKernelParameters, "yukawa");


    int N=10000;
    int verbosity=0;

    struct particles *sources = NULL;
    struct particles *targets = NULL;
    int *particleOrder = NULL;
    double *potential1 = NULL, *potential2 = NULL, *potential3 = NULL, *potential_direct = NULL;
    double potential_engy = 0;
    double potential_engy_direct = 0;

    sources = malloc(sizeof(struct particles));
    targets = malloc(sizeof(struct particles));
    potential1 = malloc(sizeof(double) * N);
    potential2 = malloc(sizeof(double) * N);
    potential3 = malloc(sizeof(double) * N);
    potential_direct = malloc(sizeof(double) * N);
    particleOrder = malloc(sizeof(int) * N);

    targets->num = 1; //single target
    targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));
    targets->order = malloc(targets->num*sizeof(int));


    sources->num = N; // 10,000 sources
    sources->x = malloc(sources->num*sizeof(double));
    sources->y = malloc(sources->num*sizeof(double));
    sources->z = malloc(sources->num*sizeof(double));
    sources->q = malloc(sources->num*sizeof(double));
    sources->w = malloc(sources->num*sizeof(double));
    sources->order = malloc(sources->num*sizeof(int));


    for (int i=0; i<targets->num; i++){
        // single target at origin
        targets->x[i]=0.0;
        targets->y[i]=0.0;
        targets->z[i]=0.0;
        targets->q[i]=1.0;
        targets->order[i] = i;
        }

    srand(1);
    for (int i=0; i<sources->num; i++){
        // 10,000 randomly distributed sources in the [-1,1] box
        sources->x[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->y[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->z[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->q[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->w[i]=((double)rand()/(double)(RAND_MAX));
        sources->order[i] = i;
        }


    int max_per_leaf=100;
    int max_per_batch=100;
    double time_tree[9];

    char *kernelName, *singularityHandling, *approximationName;


    int tree_type=1; // particle-cluster
    double kappa=0.5;
    SetKernelParameters(CoulombKernel, &kappa);
    SetKernelParameters(YukawaKernel, &kappa);

    // 3 parameter sets.  Set 2 increases order, set 3 reduces MAC.  Both should be more accurate than set 1.
    int order1=4;
    double theta1=0.8;

    int order2=5;
    double theta2=0.8;

    int order3=4;
    double theta3=0.5;

    /***********************************************/
    /******************* Test 1 ********************/
    /***********************************************/
    /********** lagrange-coulomb-skipping **********/
    /***********************************************/
    approximationName="lagrange";
    kernelName="coulomb";
    double sizeCheckFactor=0.0;
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential1[i]=0.0;
        potential2[i]=0.0;
        potential3[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, CoulombKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order1, theta1, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential1, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order2, theta2, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential2, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order3, theta3, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential3, time_tree, sizeCheckFactor, verbosity);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: lagrange-coulomb-skipping", \
                err2 < err1);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: lagrange-coulomb-skipping", \
                err3 < err1);
    }



    /***********************************************/
    /******************* Test 2 ********************/
    /***********************************************/
    /********* lagrange-coulomb-subtraction ********/
    /***********************************************/
    approximationName="lagrange";
    kernelName="coulomb";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential1[i]=0.0;
        potential2[i]=0.0;
        potential3[i]=0.0;
        potential_direct[i]=0.0;
    }



    directdriver(sources, targets, CoulombKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order1, theta1, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential1, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order2, theta2, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential2, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order3, theta3, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential3, time_tree, sizeCheckFactor, verbosity);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: lagrange-coulomb-subtraction", \
                err2 < err1);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: lagrange-coulomb-subtraction", \
                err3 < err1);
    }


    /***********************************************/
    /******************* Test 3 ********************/
    /***********************************************/
    /*********** lagrange-yukawa-skipping **********/
    /***********************************************/
    approximationName="lagrange";
    kernelName="yukawa";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential1[i]=0.0;
        potential2[i]=0.0;
        potential3[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, YukawaKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order1, theta1, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential1, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order2, theta2, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential2, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order3, theta3, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential3, time_tree, sizeCheckFactor, verbosity);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: lagrange-yukawa-skipping", \
                err2 < err1);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: lagrange-yukawa-skipping", \
                err3 < err1);
    }

    /***********************************************/
    /******************* Test 4 ********************/
    /***********************************************/
    /********* lagrange-yukawa-subtraction *********/
    /***********************************************/
    approximationName="lagrange";
    kernelName="yukawa";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential1[i]=0.0;
        potential2[i]=0.0;
        potential3[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, YukawaKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order1, theta1, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential1, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order2, theta2, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential2, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order3, theta3, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential3, time_tree, sizeCheckFactor, verbosity);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: lagrange-yukawa-subtraction", \
                err2 < err1);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: lagrange-yukawa-subtraction", \
                err3 < err1);
    }

    /***********************************************/
    /******************* Test 5 ********************/
    /***********************************************/
    /********* hermite-coulomb-skipping ************/
    /***********************************************/
    approximationName="hermite";
    sizeCheckFactor=0.0;
    kernelName="coulomb";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential1[i]=0.0;
        potential2[i]=0.0;
        potential3[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, CoulombKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order1, theta1, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential1, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order2, theta2, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential2, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order3, theta3, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential3, time_tree, sizeCheckFactor, verbosity);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: hermite-coulomb-skipping", \
                err2 < err1);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: hermite-coulomb-skipping", \
                err3 < err1);
    }


    /***********************************************/
    /******************* Test 6 ********************/
    /***********************************************/
    /******** hermite-coulomb-subtraction **********/
    /***********************************************/
    approximationName="hermite";
    kernelName="coulomb";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential1[i]=0.0;
        potential2[i]=0.0;
        potential3[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, CoulombKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order1, theta1, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential1, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order2, theta2, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential2, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order3, theta3, max_per_leaf, max_per_batch,
               CoulombKernel, singularityHandling, approximationName, tree_type,
               potential3, time_tree, sizeCheckFactor, verbosity);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: hermite-coulomb-subtraction", \
                err2 < err1);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: hermite-coulomb-subtraction", \
                err3 < err1);
    }

    /***********************************************/
    /******************* Test 7 ********************/
    /***********************************************/
    /********** hermite-yukawa-skipping ************/
    /***********************************************/
    approximationName="hermite";
    kernelName="yukawa";
    singularityHandling="skipping";
    for (int i=0; i<targets->num; i++){
        potential1[i]=0.0;
        potential2[i]=0.0;
        potential3[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, YukawaKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order1, theta1, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential1, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order2, theta2, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential2, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order3, theta3, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential3, time_tree, sizeCheckFactor, verbosity);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: hermite-yukawa-skipping", \
                err2 < err1);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: hermite-yukawa-skipping", \
                err3 < err1);
    }

    /***********************************************/
    /******************* Test 8 ********************/
    /***********************************************/
    /********* hermite-yukawa-subtraction **********/
    /***********************************************/
    approximationName="hermite";
    kernelName="yukawa";
    singularityHandling="subtraction";
    for (int i=0; i<targets->num; i++){
        potential1[i]=0.0;
        potential2[i]=0.0;
        potential3[i]=0.0;
        potential_direct[i]=0.0;
    }

    directdriver(sources, targets, YukawaKernel, singularityHandling, approximationName,
                 potential_direct, time_tree);

    treedriver(sources, targets, order1, theta1, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential1, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order2, theta2, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential2, time_tree, sizeCheckFactor, verbosity);

    treedriver(sources, targets, order3, theta3, max_per_leaf, max_per_batch,
               YukawaKernel, singularityHandling, approximationName, tree_type,
               potential3, time_tree, sizeCheckFactor, verbosity);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: hermite-yukawa-subtraction", \
                err2 < err1);
        mu_assert("TEST FAILED: increasing order didn't improve accuracy for: hermite-yukawa-subtraction", \
                err3 < err1);
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

    free(potential1);
    free(potential2);
    free(potential3);
    free(potential_direct);
    return 0;
}



// Run all the tests
static char *all_tests() {
    mu_run_test(test_direct_sum_on_10_particles);
    printf("Completed test_direct_sum_on_10_particles().\n");
    mu_run_test(test_treecode_on_100_particles);
    printf("Completed test_treecode_on_100_particles().\n");
    mu_run_test(test_treecode_on_1_target_10000_sources);
    printf("Completed test_treecode_on_1_target_10000_sources().\n");
    mu_run_test(test_treecode_parameters_on_1_target_10000_sources);
    printf("Completed test_treecode_parameters_on_1_target_10000_sources().\n");
    return 0;
}

// Run one test
static char *run_one_test(int i) {
    if (i==0){
        mu_run_test(test_direct_sum_on_10_particles);
        printf("Completed test_direct_sum_on_10_particles().\n");
    }else if(i==1){
        mu_run_test(test_treecode_on_100_particles);
        printf("Completed test_treecode_on_100_particles().\n");
    }else if(i==2){
        mu_run_test(test_treecode_on_1_target_10000_sources);
        printf("Completed test_treecode_on_1_target_10000_sources().\n");
    }else if(i==3){
        mu_run_test(test_treecode_parameters_on_1_target_10000_sources);
        printf("Completed test_treecode_parameters_on_1_target_10000_sources().\n");
    }else{
        printf("Incorrect test number.  Exiting.");
        exit(-1);
    }
    return 0;
}

int main(int argc, char **argv) {
    int rc, rank, numProcs;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

    char *result = run_one_test(atoi(argv[1]));
    if (result != 0) {
        printf("%s\n", result);
    }


    MPI_Finalize();
    return result != 0;
}
