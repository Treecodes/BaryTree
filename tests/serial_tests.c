#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <float.h>
#include <string.h>
#include <mpi.h>

#include "../src/utilities/array.h"

#include "../src/particles/struct_particles.h"

#include "../src/run_params/struct_run_params.h"
#include "../src/run_params/run_params.h"

#include "../src/interface/BaryTreeInterface.h"
#include "../src/drivers/treedriver.h"
#include "../src/drivers/directdriver.h"

#include "minunit.h"


int tests_run = 0;

static char *test_direct_sum_on_10_particles()
{
    struct RunParams *run_params = NULL;
    double time_tree[13];

    int N = 10;
    int verbosity = 1;

    struct Particles *sources = NULL;
    struct Particles *targets = NULL;
    double *potential = NULL;
    double potential_engy = 0;

    sources = malloc(sizeof(struct Particles));
    targets = malloc(sizeof(struct Particles));
    potential = malloc(sizeof(double) * N);

    targets->num = N;
    targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));

    sources->num = N;
    sources->x = malloc(sources->num*sizeof(double));
    sources->y = malloc(sources->num*sizeof(double));
    sources->z = malloc(sources->num*sizeof(double));
    sources->q = malloc(sources->num*sizeof(double));
    sources->w = malloc(sources->num*sizeof(double));


    for (int i=0; i<targets->num; i++){

        targets->x[i]=1.0*i;
        targets->y[i]=1.0*i;
        targets->z[i]=1.0*i;
        targets->q[i]=1.0*i;

        sources->x[i]=1.0*i;
        sources->y[i]=1.0*i;
        sources->z[i]=1.0*i;
        sources->q[i]=1.0*i;
        sources->w[i]=1.0;
    }

    memset(potential, 0, targets->num * sizeof(double));

    double kernelParams[1] = {0.5};
    int numKernelParams = 1;

    RunParams_Setup(&run_params,
                    COULOMB, numKernelParams, kernelParams, NO_APPROX, SKIPPING, NO_COMPUTE_TYPE,
                    0, 0, 0, 0, 0, -1, verbosity);

    directdriver(sources, targets, run_params, potential, time_tree);

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

        mu_assert("TEST FAILED: Direct sum potential not correct for coulomb kernel with skipping", \
                  fabs(potential[i] - trueValue) < 1e-10);
    }


    memset(potential, 0, targets->num * sizeof(double));

    kernelParams[0] = 5.5;
    double kappa2 = kernelParams[0] * kernelParams[0];

    RunParams_Setup(&run_params,
                    COULOMB, numKernelParams, kernelParams, NO_APPROX, SUBTRACTION, NO_COMPUTE_TYPE,
                    0, 0, 0, 0, 0, -1, verbosity);

    directdriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        double trueValue=2.0 * M_PI * kappa2 * i;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa2))/(r);
        }
        for (int j=i+1; j<targets->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa2))/(r);
        }

        mu_assert("TEST FAILED: Direct sum potential not correct for coulomb kernel with subtraction", \
                  fabs(potential[i] - trueValue) < 1e-10);
    }


    memset(potential, 0, targets->num * sizeof(double));

    kernelParams[0] = 0.5;
    double kappa = kernelParams[0];

    RunParams_Setup(&run_params,
                    YUKAWA, numKernelParams, kernelParams, NO_APPROX, SKIPPING, NO_COMPUTE_TYPE,
                    0, 0, 0, 0, 0, -1, verbosity);

    directdriver(sources, targets, run_params, potential, time_tree);

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

        mu_assert("TEST FAILED: Direct sum potential not correct for yukawa kernel with skipping", \
                  fabs(potential[i] - trueValue) < 1e-10);
    }


    memset(potential, 0, targets->num * sizeof(double));

    kernelParams[0] = 0.5;
    kappa = kernelParams[0];

    RunParams_Setup(&run_params,
                    YUKAWA, numKernelParams, kernelParams, NO_APPROX, SUBTRACTION, NO_COMPUTE_TYPE,
                    0, 0, 0, 0, 0, -1, verbosity);

    directdriver(sources, targets, run_params, potential, time_tree);

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

        mu_assert("TEST FAILED: Direct sum potential not correct for yukawa kernel with subtraction", \
                  fabs(potential[i] - trueValue) < 1e-10);
    }



    free(sources->x);
    free(sources->y);
    free(sources->z);
    free(sources->q);
    free(sources->w);
    free(sources);

    free(targets->x);
    free(targets->y);
    free(targets->z);
    free(targets->q);
    free(targets);

    free(potential);

    RunParams_Free(&run_params);

    return 0;
}


static char *test_treecode_on_100_particles()
{
    struct RunParams *run_params = NULL;
    double time_tree[13];

    int verbosity = 1;
    int N = 100;

    double beta = -1.0;

    struct Particles *sources = NULL;
    struct Particles *targets = NULL;
    double *potential = NULL;
    double potential_engy = 0;

    sources = malloc(sizeof(struct Particles));
    targets = malloc(sizeof(struct Particles));
    potential = malloc(sizeof(double) * N);

    targets->num = N;
    targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));

    sources->num = N;
    sources->x = malloc(sources->num*sizeof(double));
    sources->y = malloc(sources->num*sizeof(double));
    sources->z = malloc(sources->num*sizeof(double));
    sources->q = malloc(sources->num*sizeof(double));
    sources->w = malloc(sources->num*sizeof(double));


    for (int i=0; i<targets->num; i++){
        targets->x[i]=1.0*i;
        targets->y[i]=1.0*i;
        targets->z[i]=1.0*i;
        targets->q[i]=1.0*i;

        sources->x[i]=1.0*i;
        sources->y[i]=1.0*i;
        sources->z[i]=1.0*i;
        sources->q[i]=1.0*i;
        sources->w[i]=1.0;
    }


    int max_per_source_leaf = 3;
    int max_per_target_leaf = 3;

    int degree = 2;
    double theta = 0.7;
    double size_check = 1.0;

    int num_kernel_params = 1;
    double kernel_params[1] = {0.5};
    double kappa = kernel_params[0];

    RunParams_Setup(&run_params,
                    NO_KERNEL, num_kernel_params, kernel_params, NO_APPROX, NO_SINGULARITY, PARTICLE_CLUSTER,
                    theta, degree, max_per_source_leaf, max_per_target_leaf, size_check, beta, verbosity);


    /***********************************************/
    /******************* Test 1 ********************/
    /***********************************************/
    /********** lagrange-coulomb-skipping **********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SKIPPING;
    run_params->approximation = LAGRANGE;

    treedriver(sources, targets, run_params, potential, time_tree);

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
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-coulomb-skipping", \
                  fabs(potential[i] - trueValue)/fabs(trueValue) < 3e-3);
    }


    /***********************************************/
    /******************* Test 2 ********************/
    /***********************************************/
    /********* lagrange-coulomb-subtraction ********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = LAGRANGE;

    treedriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        double trueValue=2.0 * M_PI * kappa * kappa * i;
        for (int j=0; j<i; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa*kappa))/(r);
        }
        for (int j=i+1; j<targets->num; j++){
            double r = abs(j-i)*sqrt(3);
            trueValue += (j - i*exp(-r*r/kappa*kappa))/(r);
        }

        if (verbosity>0) printf("trueValue = %1.8e\n", trueValue);
        if (verbosity>0) printf("computedValue = %1.8e\n", potential[i]);
        if (verbosity>0) printf("relative error = %1.8e\n\n", fabs(potential[i] - trueValue)/fabs(trueValue));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-coulomb-subtraction", \
                  fabs(potential[i] - trueValue)/fabs(trueValue) < 2e-2);
    }


    /***********************************************/
    /******************* Test 3 ********************/
    /***********************************************/
    /*********** lagrange-yukawa-skipping **********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SKIPPING; 
    run_params->approximation = LAGRANGE;

    treedriver(sources, targets, run_params, potential, time_tree);

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
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-yukawa-skipping", \
                  fabs(potential[i] - trueValue)/fabs(trueValue) < 8e-3);
    }


    /***********************************************/
    /******************* Test 4 ********************/
    /***********************************************/
    /********* lagrange-yukawa-subtraction *********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = LAGRANGE;

    treedriver(sources, targets, run_params, potential, time_tree);

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
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-yukawa-subtraction", \
                  fabs(potential[i] - trueValue) < 2e-2);
    }


    /***********************************************/
    /******************* Test 5 ********************/
    /***********************************************/
    /********* hermite-coulomb-skipping ************/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SKIPPING;
    run_params->approximation = HERMITE;

    treedriver(sources, targets, run_params, potential, time_tree);

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
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-coulomb-skipping", \
                  fabs(potential[i] - trueValue)/fabs(trueValue) < 3e-4);
    }


    /***********************************************/
    /******************* Test 6 ********************/
    /***********************************************/
    /******** hermite-coulomb-subtraction **********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = HERMITE;

    treedriver(sources, targets, run_params, potential, time_tree);

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
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-coulomb-subtraction", \
                  fabs(potential[i] - trueValue)/fabs(trueValue) < 2e-2);
    }


    /***********************************************/
    /******************* Test 7 ********************/
    /***********************************************/
    /********** hermite-yukawa-skipping ************/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SKIPPING;
    run_params->approximation = HERMITE;

    treedriver(sources, targets, run_params, potential, time_tree);

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
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-yukawa-skipping", \
                  fabs(potential[i] - trueValue)/fabs(trueValue) < 5e-4);
    }


    /***********************************************/
    /******************* Test 8 ********************/
    /***********************************************/
    /********* hermite-yukawa-subtraction **********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = HERMITE;

    treedriver(sources, targets, run_params, potential, time_tree);

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
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-yukawa-subtraction", \
                  fabs(potential[i] - trueValue) < 2e-3);

    }

    free(sources->x);
    free(sources->y);
    free(sources->z);
    free(sources->q);
    free(sources->w);
    free(sources);

    free(targets->x);
    free(targets->y);
    free(targets->z);
    free(targets->q);
    free(targets);

    free(potential);

    RunParams_Free(&run_params);

    return 0;
}


static char *test_treecode_on_1_target_10000_sources()
{
    struct RunParams *run_params = NULL;
    double time_tree[13];

    int verbosity = 1;
    int N = 10000;

    double beta = -1.0;

    struct Particles *sources = NULL;
    struct Particles *targets = NULL;
    double *potential = NULL, *potential_direct = NULL;
    double potential_engy = 0;
    double potential_engy_direct = 0;

    sources = malloc(sizeof(struct Particles));
    targets = malloc(sizeof(struct Particles));
    potential = malloc(sizeof(double) * N);
    potential_direct = malloc(sizeof(double) * N);

    targets->num = 1; //single target
    targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));


    sources->num = N; // 10,000 sources
    sources->x = malloc(sources->num*sizeof(double));
    sources->y = malloc(sources->num*sizeof(double));
    sources->z = malloc(sources->num*sizeof(double));
    sources->q = malloc(sources->num*sizeof(double));
    sources->w = malloc(sources->num*sizeof(double));


    for (int i=0; i<targets->num; i++){
        // single target at origin
        targets->x[i]=0.0;
        targets->y[i]=0.0;
        targets->z[i]=0.0;
        targets->q[i]=1.0;
    }

    srand(1);
    for (int i=0; i<sources->num; i++){
        // 10,000 randomly distributed sources in the [-1,1] box
        sources->x[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->y[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->z[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->q[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->w[i]=((double)rand()/(double)(RAND_MAX));
    }


    int max_per_source_leaf = 100;
    int max_per_target_leaf = 100;

    int degree = 4;
    double theta = 0.8;
    double size_check = 0.0;

    int num_kernel_params = 1;
    double kernel_params[1] = {0.5};

    RunParams_Setup(&run_params,
                    NO_KERNEL, num_kernel_params, kernel_params, NO_APPROX, NO_SINGULARITY, PARTICLE_CLUSTER,
                    theta, degree, max_per_source_leaf, max_per_target_leaf, size_check, beta, verbosity);


    /***********************************************/
    /******************* Test 1 ********************/
    /***********************************************/
    /********** lagrange-coulomb-skipping **********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SKIPPING;
    run_params->approximation = LAGRANGE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);
    treedriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nlagrange-coulomb-skipping\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])
                                                         / fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-coulomb-skipping", \
                  fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-4);
    }


    /***********************************************/
    /******************* Test 2 ********************/
    /***********************************************/
    /********* lagrange-coulomb-subtraction ********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = LAGRANGE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);
    treedriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nlagrange-coulomb-subtraction\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])
                                                         / fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-coulomb-subtraction", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-4);
    }


    /***********************************************/
    /******************* Test 3 ********************/
    /***********************************************/
    /*********** lagrange-yukawa-skipping **********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SKIPPING;
    run_params->approximation = LAGRANGE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);
    treedriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nlagrange-yukawa-skipping\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])
                                                         / fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-yukawa-skipping", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-4);
    }


    /***********************************************/
    /******************* Test 4 ********************/
    /***********************************************/
    /********* lagrange-yukawa-subtraction *********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = LAGRANGE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);
    treedriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nlagrange-yukawa-subtraction\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])
                                                         / fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: lagrange-yukawa-subtraction", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-4);
    }


    /***********************************************/
    /******************* Test 5 ********************/
    /***********************************************/
    /********* hermite-coulomb-skipping ************/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SKIPPING;
    run_params->approximation = HERMITE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);
    treedriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nhermite-coulomb-skipping\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])
                                                         / fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-coulomb-skipping", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-5);
    }


    /***********************************************/
    /******************* Test 6 ********************/
    /***********************************************/
    /******** hermite-coulomb-subtraction **********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = HERMITE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);
    treedriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nhermite-coulomb-subtraction\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])
                                                         / fabs(potential_direct[i]));
        mu_assert("TEST FAILED: Treecode potential not correct for: hermite-coulomb-subtraction", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-5);
    }


    /***********************************************/
    /******************* Test 7 ********************/
    /***********************************************/
    /********** hermite-yukawa-skipping ************/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SKIPPING;
    run_params->approximation = HERMITE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);
    treedriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nhermite-yukawa-skipping\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])
                                                         / fabs(potential_direct[i]));
        mu_assert("TEST FAILED: TEST FAILED: Treecode potential not correct for: hermite-yukawa-skipping", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-5);
    }


    /***********************************************/
    /******************* Test 8 ********************/
    /***********************************************/
    /********* hermite-yukawa-subtraction **********/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = HERMITE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);
    treedriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nhermite-yukawa-subtraction\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_direct[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_direct[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_direct[i])
                                                         / fabs(potential_direct[i]));
        mu_assert("TEST FAILED: TEST FAILED: Treecode potential not correct for: hermite-yukawa-subtraction", \
                fabs(potential[i] - potential_direct[i])/fabs(potential_direct[i]) < 2e-5);
    }

    free(sources->x);
    free(sources->y);
    free(sources->z);
    free(sources->q);
    free(sources->w);
    free(sources);

    free(targets->x);
    free(targets->y);
    free(targets->z);
    free(targets->q);
    free(targets);

    free(potential);
    free(potential_direct);

    RunParams_Free(&run_params);

    return 0;
}


static char *test_treecode_wrapper()
{
    struct RunParams *run_params = NULL;
    double time_tree[13];

    int verbosity = 1;
    int N = 10000;

    struct Particles *sources = NULL;
    struct Particles *targets = NULL;
    double *potential = NULL, *potential_wrapper = NULL;
    double potential_engy = 0;
    double potential_engy_direct = 0;

    sources = malloc(sizeof(struct Particles));
    targets = malloc(sizeof(struct Particles));
    potential = malloc(sizeof(double) * N);
    potential_wrapper = malloc(sizeof(double) * N);

    targets->num = 1; //single target
    targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));

    sources->num = N; // 10,000 sources
    sources->x = malloc(sources->num*sizeof(double));
    sources->y = malloc(sources->num*sizeof(double));
    sources->z = malloc(sources->num*sizeof(double));
    sources->q = malloc(sources->num*sizeof(double));
    sources->w = malloc(sources->num*sizeof(double));


    for (int i=0; i<targets->num; i++){
        // single target at origin
        targets->x[i]=0.0;
        targets->y[i]=0.0;
        targets->z[i]=0.0;
        targets->q[i]=1.0;
    }

    srand(1);
    for (int i=0; i<sources->num; i++){
        // 10,000 randomly distributed sources in the [-1,1] box
        sources->x[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->y[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->z[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->q[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->w[i]=((double)rand()/(double)(RAND_MAX));
    }


    int max_per_source_leaf = 100;
    int max_per_target_leaf = 100;

    int degree = 4;
    double theta = 0.8;
    double beta = -1.0;
    double size_check = 1.0;

    int num_kernel_params = 1;
    double kernel_params[1] = {0.5};

    RunParams_Setup(&run_params,
                    NO_KERNEL, num_kernel_params, kernel_params, NO_APPROX, NO_SINGULARITY, PARTICLE_CLUSTER,
                    theta, degree, max_per_source_leaf, max_per_target_leaf, size_check, beta, verbosity);


    /***********************************************/
    /******************* Test **********************/
    /***********************************************/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_wrapper, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SKIPPING;
    run_params->approximation = LAGRANGE;

    BaryTreeInterface(targets->num, sources->num,
                      targets->x,targets->y,targets->z,targets->q,
                      sources->x,sources->y,sources->z,sources->q,sources->w,
                      potential_wrapper, COULOMB, num_kernel_params, kernel_params,
                      SKIPPING, LAGRANGE, PARTICLE_CLUSTER,
                      theta, degree, max_per_source_leaf, max_per_target_leaf,
                      size_check, beta, verbosity);

    treedriver(sources, targets, run_params, potential, time_tree);

    for (int i=0; i<targets->num; i++){
        if (verbosity>0) printf("\nlagrange-coulomb-skipping\n");
        if (verbosity>0) printf("direct: %1.8e\n", potential_wrapper[i]);
        if (verbosity>0) printf("approx: %1.8e\n", potential[i]);
        if (verbosity>0) printf("absolute error: %1.2e\n", fabs(potential[i] - potential_wrapper[i]));
        if (verbosity>0) printf("relative error: %1.2e\n", fabs(potential[i] - potential_wrapper[i])
                                                         / fabs(potential_wrapper[i]));
        mu_assert("TEST FAILED: Treecode wrapper didn't give same results as directly calling treedriver", \
                fabs(potential[i] - potential_wrapper[i])/fabs(potential[i]) < 1e-12);
    }

    free(sources->x);
    free(sources->y);
    free(sources->z);
    free(sources->q);
    free(sources->w);
    free(sources);

    free(targets->x);
    free(targets->y);
    free(targets->z);
    free(targets->q);
    free(targets);

    free(potential);
    free(potential_wrapper);

    RunParams_Free(&run_params);

    return 0;
}


static char *test_treecode_parameters_on_1_target_5000_sources()
{
    struct RunParams *run_params = NULL;
    double time_tree[9];

    int verbosity = 1;
    int N = 5000;
    double beta = -1.0;

    struct Particles *sources = NULL;
    struct Particles *targets = NULL;
    double *potential1 = NULL, *potential2 = NULL, *potential3 = NULL, *potential_direct = NULL;
    double potential_engy = 0;
    double potential_engy_direct = 0;


    sources = malloc(sizeof(struct Particles));
    targets = malloc(sizeof(struct Particles));
    potential1 = malloc(sizeof(double) * N);
    potential2 = malloc(sizeof(double) * N);
    potential3 = malloc(sizeof(double) * N);
    potential_direct = malloc(sizeof(double) * N);

    targets->num = 1; //single target
    targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));

    sources->num = N; // 10,000 sources
    sources->x = malloc(sources->num*sizeof(double));
    sources->y = malloc(sources->num*sizeof(double));
    sources->z = malloc(sources->num*sizeof(double));
    sources->q = malloc(sources->num*sizeof(double));
    sources->w = malloc(sources->num*sizeof(double));


    for (int i=0; i<targets->num; i++){
        // single target at origin
        targets->x[i]=0.0;
        targets->y[i]=0.0;
        targets->z[i]=0.0;
        targets->q[i]=1.0;
    }

    srand(1);
    for (int i=0; i<sources->num; i++){
        // 10,000 randomly distributed sources in the [-1,1] box
        sources->x[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->y[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->z[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->q[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->w[i]=((double)rand()/(double)(RAND_MAX));
    }


    int max_per_source_leaf = 50;
    int max_per_target_leaf = 5;
    double size_check = 0.0;

    int num_kernel_params = 1;
    double kernel_params[1] = {0.5};


    RunParams_Setup(&run_params,
                    NO_KERNEL, num_kernel_params, kernel_params, NO_APPROX, NO_SINGULARITY, PARTICLE_CLUSTER,
                    0, 0, max_per_source_leaf, max_per_target_leaf, size_check, beta, verbosity);


    // 3 parameter sets.  Set 2 increases degree, set 3 reduces MAC.  Both should be more accurate than set 1.
    int degree1=3;
    double theta1=0.9;

    int degree2=6;
    double theta2=0.9;

    int degree3=3;
    double theta3=0.4;

    /***********************************************/
    /******************* Test 1 ********************/
    /***********************************************/
    /********** lagrange-coulomb-skipping **********/
    /***********************************************/
    memset(potential1, 0, targets->num * sizeof(double));
    memset(potential2, 0, targets->num * sizeof(double));
    memset(potential3, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SKIPPING;
    run_params->approximation = LAGRANGE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);

    run_params->interp_degree = degree1;
    run_params->theta = theta1;
    treedriver(sources, targets, run_params, potential1, time_tree);

    run_params->interp_degree = degree2;
    run_params->theta = theta2;
    treedriver(sources, targets, run_params, potential2, time_tree);

    run_params->interp_degree = degree3;
    run_params->theta = theta3;
    treedriver(sources, targets, run_params, potential3, time_tree);


    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing degree didn't improve accuracy for: lagrange-coulomb-skipping", \
                err2 < err1);
        mu_assert("TEST FAILED: decreasing theta didn't improve accuracy for: lagrange-coulomb-skipping", \
                err3 < err1);
    }


    /***********************************************/
    /******************* Test 2 ********************/
    /***********************************************/
    /********* lagrange-coulomb-subtraction ********/
    /***********************************************/
    memset(potential1, 0, targets->num * sizeof(double));
    memset(potential2, 0, targets->num * sizeof(double));
    memset(potential3, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = LAGRANGE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);

    run_params->interp_degree = degree1;
    run_params->theta = theta1;
    treedriver(sources, targets, run_params, potential1, time_tree);

    run_params->interp_degree = degree2;
    run_params->theta = theta2;
    treedriver(sources, targets, run_params, potential2, time_tree);

    run_params->interp_degree = degree3;
    run_params->theta = theta3;
    treedriver(sources, targets, run_params, potential3, time_tree);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing degree didn't improve accuracy for: lagrange-coulomb-subtraction", \
                err2 < err1);
        mu_assert("TEST FAILED: decreasing theta didn't improve accuracy for: lagrange-coulomb-subtraction", \
                err3 < err1);
    }


    /***********************************************/
    /******************* Test 3 ********************/
    /***********************************************/
    /*********** lagrange-yukawa-skipping **********/
    /***********************************************/
    memset(potential1, 0, targets->num * sizeof(double));
    memset(potential2, 0, targets->num * sizeof(double));
    memset(potential3, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SKIPPING;
    run_params->approximation = LAGRANGE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);

    run_params->interp_degree = degree1;
    run_params->theta = theta1;
    treedriver(sources, targets, run_params, potential1, time_tree);

    run_params->interp_degree = degree2;
    run_params->theta = theta2;
    treedriver(sources, targets, run_params, potential2, time_tree);

    run_params->interp_degree = degree3;
    run_params->theta = theta3;
    treedriver(sources, targets, run_params, potential3, time_tree);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing degree didn't improve accuracy for: lagrange-yukawa-skipping", \
                err2 < err1);
        mu_assert("TEST FAILED: decreasing theta didn't improve accuracy for: lagrange-yukawa-skipping", \
                err3 < err1);
    }


    /***********************************************/
    /******************* Test 4 ********************/
    /***********************************************/
    /********* lagrange-yukawa-subtraction *********/
    /***********************************************/
    memset(potential1, 0, targets->num * sizeof(double));
    memset(potential2, 0, targets->num * sizeof(double));
    memset(potential3, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = LAGRANGE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);

    run_params->interp_degree = degree1;
    run_params->theta = theta1;
    treedriver(sources, targets, run_params, potential1, time_tree);

    run_params->interp_degree = degree2;
    run_params->theta = theta2;
    treedriver(sources, targets, run_params, potential2, time_tree);

    run_params->interp_degree = degree3;
    run_params->theta = theta3;
    treedriver(sources, targets, run_params, potential3, time_tree);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing degree didn't improve accuracy for: lagrange-yukawa-subtraction", \
                err2 < err1);
        mu_assert("TEST FAILED: decreasing theta didn't improve accuracy for: lagrange-yukawa-subtraction", \
                err3 < err1);
    }


    /***********************************************/
    /******************* Test 5 ********************/
    /***********************************************/
    /********* hermite-coulomb-skipping ************/
    /***********************************************/
    memset(potential1, 0, targets->num * sizeof(double));
    memset(potential2, 0, targets->num * sizeof(double));
    memset(potential3, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SKIPPING;
    run_params->approximation = HERMITE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);

    run_params->interp_degree = degree1;
    run_params->theta = theta1;
    treedriver(sources, targets, run_params, potential1, time_tree);

    run_params->interp_degree = degree2;
    run_params->theta = theta2;
    treedriver(sources, targets, run_params, potential2, time_tree);

    run_params->interp_degree = degree3;
    run_params->theta = theta3;
    treedriver(sources, targets, run_params, potential3, time_tree);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        verbosity=1;
        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        verbosity=0;
        mu_assert("TEST FAILED: increasing degree didn't improve accuracy for: hermite-coulomb-skipping", \
                err2 < err1);
        mu_assert("TEST FAILED: decreasing theta didn't improve accuracy for: hermite-coulomb-skipping", \
                err3 < err1);
    }


    /***********************************************/
    /******************* Test 6 ********************/
    /***********************************************/
    /******** hermite-coulomb-subtraction **********/
    /***********************************************/
    memset(potential1, 0, targets->num * sizeof(double));
    memset(potential2, 0, targets->num * sizeof(double));
    memset(potential3, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = HERMITE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);

    run_params->interp_degree = degree1;
    run_params->theta = theta1;
    treedriver(sources, targets, run_params, potential1, time_tree);

    run_params->interp_degree = degree2;
    run_params->theta = theta2;
    treedriver(sources, targets, run_params, potential2, time_tree);

    run_params->interp_degree = degree3;
    run_params->theta = theta3;
    treedriver(sources, targets, run_params, potential3, time_tree);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing degree didn't improve accuracy for: hermite-coulomb-subtraction", \
                err2 < err1);
        mu_assert("TEST FAILED: decreasing theta didn't improve accuracy for: hermite-coulomb-subtraction", \
                err3 < err1);
    }


    /***********************************************/
    /******************* Test 7 ********************/
    /***********************************************/
    /********** hermite-yukawa-skipping ************/
    /***********************************************/
    memset(potential1, 0, targets->num * sizeof(double));
    memset(potential2, 0, targets->num * sizeof(double));
    memset(potential3, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SKIPPING;
    run_params->approximation = HERMITE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);

    run_params->interp_degree = degree1;
    run_params->theta = theta1;
    treedriver(sources, targets, run_params, potential1, time_tree);

    run_params->interp_degree = degree2;
    run_params->theta = theta2;
    treedriver(sources, targets, run_params, potential2, time_tree);

    run_params->interp_degree = degree3;
    run_params->theta = theta3;
    treedriver(sources, targets, run_params, potential3, time_tree);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing degree didn't improve accuracy for: hermite-yukawa-skipping", \
                err2 < err1);
        mu_assert("TEST FAILED: decreasing theta didn't improve accuracy for: hermite-yukawa-skipping", \
                err3 < err1);
    }


    /***********************************************/
    /******************* Test 8 ********************/
    /***********************************************/
    /********* hermite-yukawa-subtraction **********/
    /***********************************************/
    memset(potential1, 0, targets->num * sizeof(double));
    memset(potential2, 0, targets->num * sizeof(double));
    memset(potential3, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = YUKAWA;
    run_params->singularity   = SUBTRACTION;
    run_params->approximation = HERMITE;

    directdriver(sources, targets, run_params, potential_direct, time_tree);

    run_params->interp_degree = degree1;
    run_params->theta = theta1;
    treedriver(sources, targets, run_params, potential1, time_tree);

    run_params->interp_degree = degree2;
    run_params->theta = theta2;
    treedriver(sources, targets, run_params, potential2, time_tree);

    run_params->interp_degree = degree3;
    run_params->theta = theta3;
    treedriver(sources, targets, run_params, potential3, time_tree);

    for (int i=0; i<targets->num; i++){
        double err1 = fabs(potential1[i] - potential_direct[i]);
        double err2 = fabs(potential2[i] - potential_direct[i]);
        double err3 = fabs(potential3[i] - potential_direct[i]);

        if (verbosity>0) printf("err1 = %1.4e\n", err1);
        if (verbosity>0) printf("err2 = %1.4e\n", err2);
        if (verbosity>0) printf("err3 = %1.4e\n", err3);
        mu_assert("TEST FAILED: increasing degree didn't improve accuracy for: hermite-yukawa-subtraction", \
                err2 < err1);
        mu_assert("TEST FAILED: decreasing theta didn't improve accuracy for: hermite-yukawa-subtraction", \
                err3 < err1);
    }

    free(sources->x);
    free(sources->y);
    free(sources->z);
    free(sources->q);
    free(sources->w);
    free(sources);

    free(targets->x);
    free(targets->y);
    free(targets->z);
    free(targets->q);
    free(targets);

    free(potential1);
    free(potential2);
    free(potential3);
    free(potential_direct);

    RunParams_Free(&run_params);

    return 0;
}





static char *test_BLDTT()
{
    struct RunParams *run_params = NULL;
    double time_tree[13];

    int verbosity = 1;
    int N = 5000;

    struct Particles *sources = NULL;
    struct Particles *targets = NULL;
    double *potential = NULL, *potential_direct = NULL;
    double potential_engy = 0;
    double potential_engy_direct = 0;

    sources = malloc(sizeof(struct Particles));
    targets = malloc(sizeof(struct Particles));
    potential = malloc(sizeof(double) * N);
    potential_direct = malloc(sizeof(double) * N);

    targets->num = N;
    targets->x = malloc(targets->num*sizeof(double));
    targets->y = malloc(targets->num*sizeof(double));
    targets->z = malloc(targets->num*sizeof(double));
    targets->q = malloc(targets->num*sizeof(double));

    sources->num = N;
    sources->x = malloc(sources->num*sizeof(double));
    sources->y = malloc(sources->num*sizeof(double));
    sources->z = malloc(sources->num*sizeof(double));
    sources->q = malloc(sources->num*sizeof(double));
    sources->w = malloc(sources->num*sizeof(double));


    srand(1);
    for (int i=0; i<sources->num; i++){
        // 10,000 randomly distributed sources in the [-1,1] box
        targets->x[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        targets->y[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        targets->z[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        targets->q[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;

        sources->x[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->y[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->z[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->q[i]=((double)rand()/(double)(RAND_MAX)) * 2. - 1.;
        sources->w[i]=((double)rand()/(double)(RAND_MAX));
    }


    int max_per_source_leaf = 20;
    int max_per_target_leaf = 20;

    int degree = 3;
    double theta = 0.9;
    double beta = -1.0;
    double size_check = 1.0;

    int num_kernel_params = 1;
    double kernel_params[1] = {0.5};

    RunParams_Setup(&run_params,
                    NO_KERNEL, num_kernel_params, kernel_params, NO_APPROX, NO_SINGULARITY, CLUSTER_CLUSTER,
                    theta, degree, max_per_source_leaf, max_per_target_leaf, size_check, beta, verbosity);


    /***********************************************/
    /******************* Test **********************/
    /***********************************************/
    /***********************************************/
    memset(potential, 0, targets->num * sizeof(double));
    memset(potential_direct, 0, targets->num * sizeof(double));

    run_params->kernel        = COULOMB;
    run_params->singularity   = SKIPPING;
    run_params->approximation = LAGRANGE;


    directdriver(sources, targets, run_params, potential_direct, time_tree);

    treedriver(sources, targets, run_params, potential, time_tree);

    double cumulative_potential_bldtt=0.0;
    double cumulative_potential=0.0;
    double error;

    for (int i=0; i<targets->num; i++){

        cumulative_potential_bldtt += potential[i];
        cumulative_potential += potential_direct[i];

    }

    error=fabs(cumulative_potential_bldtt - cumulative_potential)/fabs(cumulative_potential);

    if (verbosity>-1) printf("direct: %1.8e\n", cumulative_potential);
    if (verbosity>-1) printf("approx: %1.8e\n", cumulative_potential_bldtt);
    if (verbosity>-1) printf("rel. err.: %1.8e\n", error);

    mu_assert("TEST FAILED: Cluster-cluster didn't give same results as direct", \
            error < 6e-3);

    free(sources->x);
    free(sources->y);
    free(sources->z);
    free(sources->q);
    free(sources->w);
    free(sources);

    free(targets->x);
    free(targets->y);
    free(targets->z);
    free(targets->q);
    free(targets);

    free(potential);
    free(potential_direct);

    RunParams_Free(&run_params);

    return 0;
}



// Run all the tests
static char *all_tests()
{
    mu_run_test(test_direct_sum_on_10_particles);
    printf("Completed test_direct_sum_on_10_particles().\n");
    mu_run_test(test_treecode_on_100_particles);
    printf("Completed test_treecode_on_100_particles().\n");
    mu_run_test(test_treecode_on_1_target_10000_sources);
    printf("Completed test_treecode_on_1_target_10000_sources().\n");
    mu_run_test(test_treecode_parameters_on_1_target_5000_sources);
    printf("Completed test_treecode_parameters_on_1_target_5000_sources().\n");
    return 0;
}

// Run one test
static char *run_one_test(int i)
{
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
        mu_run_test(test_treecode_parameters_on_1_target_5000_sources);
        printf("Completed test_treecode_parameters_on_1_target_5000_sources().\n");
    }else if (i==4){
        mu_run_test(test_treecode_wrapper);
        printf("Completed test_treecode_wrapper().\n");
    }else if (i==5){
        mu_run_test(test_BLDTT);
        printf("Completed test_BLDTT().\n");
    }else{
        printf("Incorrect test number.  Exiting.\n");
        exit(1);
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
