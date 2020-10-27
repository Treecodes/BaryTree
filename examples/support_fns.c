#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "../src/utilities/tools.h"
#include "../src/run_params/run_params.h"
#include "../src/run_params/struct_run_params.h"

#include "support_fns.h"

static double erfinv (double x);


void Params_Parse(FILE *fp, struct RunParams **run_params, int *N, int *M, int *run_direct, int *slice,
                  double *xyz_limits, DISTRIBUTION *distribution, PARTITION *partition)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    /* BaryTree params */
    int verbosity = 0;
    int interp_degree = 5;
    double theta = 0.5; 
    double beta = -1.0;
    int max_per_source_leaf = 500;
    int max_per_target_leaf = 500; 
    double size_check_factor = 1.0;
    
    char kernel_string[256]        = "COULOMB";
    char singularity_string[256]   = "SKIPPING";
    char approximation_string[256] = "LAGRANGE";
    char compute_type_string[256]  = "PARTICLE_CLUSTER";
    char run_direct_string[256]    = "OFF";
    char distribution_string[256]  = "UNIFORM";
    char partition_string[256]     = "RCB";

    KERNEL kernel;
    SINGULARITY singularity;
    APPROXIMATION approximation;
    COMPUTE_TYPE compute_type;

    int num_kernel_params = 0;
    double kernel_params[32];

    /* random_cube_example params */
    *N = 10000; 
    *M = 10000;
    *slice = 1;
    
    xyz_limits[0] = -1.;
    xyz_limits[1] =  1.;
    xyz_limits[2] = -1.;
    xyz_limits[3] =  1.;
    xyz_limits[4] = -1.;
    xyz_limits[5] =  1.;


    char c[256], c1[256], c2[256];

    while (fgets(c, 256, fp) != NULL) {
        sscanf(c, "%s %s", c1, c2);
    
        /* Parameters for the RunParam struct */
        if (strcmp(c1, "degree") == 0) {
            interp_degree = atoi(c2);

        } else if (strcmp(c1, "theta") == 0) {
            theta = atof(c2);

        } else if (strcmp(c1, "beta") == 0) {
            beta = atof(c2);

        } else if (strcmp(c1, "max_per_source_leaf") == 0) {
            max_per_source_leaf = atoi(c2);

        } else if (strcmp(c1, "max_per_target_leaf") == 0) {
            max_per_target_leaf = atoi(c2);

        } else if (strcmp(c1, "kernel_name") == 0) {
            strcpy(kernel_string, c2);
        
        } else if (strcmp(c1, "kernel_params") == 0) {
            char *k_param = strtok(c2, " ,");
            while (k_param != NULL) {
                kernel_params[num_kernel_params] = atof(k_param);
                num_kernel_params++;
                k_param = strtok(NULL, " ,");
            }
        
        } else if (strcmp(c1, "singularity") == 0) {
            strcpy(singularity_string, c2);
        
        } else if (strcmp(c1, "approximation") == 0) {
            strcpy(approximation_string, c2);
        
        } else if (strcmp(c1, "compute_type") == 0) {
            strcpy(compute_type_string, c2);
        
        } else if (strcmp(c1, "size_check") == 0) {
            size_check_factor = atof(c2);
        
        } else if (strcmp(c1, "verbosity") == 0) {
            verbosity = atoi(c2);

        /* Other run parameters */
        } else if (strcmp(c1, "num_particles") == 0) {
            *N = atoi(c2);
            *M = atoi(c2);

        } else if (strcmp(c1, "num_sources") == 0) {
            *N = atoi(c2);

        } else if (strcmp(c1, "num_targets") == 0) {
            *M = atoi(c2);

        } else if (strcmp(c1, "run_direct") == 0) {
            strcpy(run_direct_string, c2);

        } else if (strcmp(c1, "slice") == 0) {
            *slice = atoi(c2);
            
        } else if (strcmp(c1, "distribution") == 0) {
            strcpy(distribution_string, c2);

        } else if (strcmp(c1, "partition") == 0) {
            strcpy(partition_string, c2);

        } else {
            if (rank == 0) {
                printf("[random cube example] ERROR! Undefined token \"%s\". Exiting.\n", c1);
            }
            exit(1);
        }
    }
    

    /* Validating tokens for RunParam struct */
    if (strcasecmp(kernel_string, "COULOMB") == 0) {
        kernel = COULOMB;

    } else if (strcasecmp(kernel_string, "YUKAWA") == 0) {
        kernel = YUKAWA;
    
    } else if ((strcasecmp(kernel_string, "REGULARIZED_COULOMB") == 0)
            || (strcasecmp(kernel_string, "REGULARIZED-COULOMB") == 0)) {
        kernel = REGULARIZED_COULOMB;
    
    } else if ((strcasecmp(kernel_string, "REGULARIZED_YUKAWA") == 0)
            || (strcasecmp(kernel_string, "REGULARIZED-YUKAWA") == 0)) {
        kernel = REGULARIZED_YUKAWA;
    
    } else if (strcasecmp(kernel_string, "ATAN") == 0) {
        kernel = ATAN;

    } else if (strcasecmp(kernel_string, "MQ") == 0) {
        kernel = MQ;
    
    } else if (strcasecmp(kernel_string, "TCF") == 0) {
        kernel = TCF;
    
    } else if (strcasecmp(kernel_string, "USER") == 0) {
        kernel = USER;

    } else if (strcasecmp(kernel_string, "DCF") == 0) {
        kernel = DCF;

    } else if ((strcasecmp(kernel_string, "SIN_OVER_R") == 0)
            || (strcasecmp(kernel_string, "SIN-OVER-R") == 0)) {
        kernel = SIN_OVER_R;
    
    } else {
        if (rank == 0) {
            printf("[random cube example] ERROR! Undefined kernel token \"%s\". Exiting.\n",
                   kernel_string);
        }
        exit(1);
    }


    if (strcasecmp(singularity_string, "SKIPPING") == 0) {
        singularity = SKIPPING;

    } else if (strcasecmp(singularity_string, "SUBTRACTION") == 0) {
        singularity = SUBTRACTION;

    } else {
        if (rank == 0) {
            printf("[random cube example] ERROR! Undefined singularity token \"%s\". Exiting.\n",
                   singularity_string);
        }
        exit(1);
    }


    if (strcasecmp(approximation_string, "LAGRANGE") == 0) {
        approximation = LAGRANGE;

    } else if (strcasecmp(approximation_string, "HERMITE") == 0) {
        approximation = HERMITE;
    
    } else {
        if (rank == 0) {
            printf("[random cube example] ERROR! Undefined approximation token \"%s\". Exiting.\n",
                   approximation_string);
        }
        exit(1);
    }


    if ((strcasecmp(compute_type_string, "CLUSTER_PARTICLE") == 0)
     || (strcasecmp(compute_type_string, "CLUSTER-PARTICLE") == 0)) {
        compute_type = CLUSTER_PARTICLE;

    } else if ((strcasecmp(compute_type_string, "PARTICLE_CLUSTER") == 0)
            || (strcasecmp(compute_type_string, "PARTICLE-CLUSTER") == 0)) {
        compute_type = PARTICLE_CLUSTER;

    } else if ((strcasecmp(compute_type_string, "CLUSTER_CLUSTER") == 0)
            || (strcasecmp(compute_type_string, "CLUSTER-CLUSTER") == 0)) {
        compute_type = CLUSTER_CLUSTER;

    } else {
        if (rank == 0) {
            printf("[random cube example] ERROR! Undefined compute_type token \"%s\". Exiting.\n",
                   compute_type_string);
        }
        exit(1);
    }
    
    
    /* Validating other tokens */
    if ((strcasecmp(run_direct_string, "ON")  == 0)
     || (strcasecmp(run_direct_string, "YES") == 0)
     || (strcasecmp(run_direct_string, "1")   == 0)) {
        *run_direct = 1;
        
    } else if ((strcasecmp(run_direct_string, "OFF") == 0)
            || (strcasecmp(run_direct_string, "NO")  == 0)
            || (strcasecmp(run_direct_string, "0")   == 0)) {
        *run_direct = 0;
    } else {
        if (rank == 0) {
            printf("[random cube example] ERROR! Undefined run direct token \"%s\". Exiting.\n",
                   run_direct_string);
        }
        exit(1);
    }
    
    
    if (strcasecmp(distribution_string, "UNIFORM") == 0) {
        *distribution = UNIFORM;
        
    } else if ((strcasecmp(distribution_string, "GAUSSIAN") == 0)
            || (strcasecmp(distribution_string, "NORMAL") == 0)) {
        *distribution = GAUSSIAN;
        
    } else if (strcasecmp(distribution_string, "EXPONENTIAL") == 0) {
        *distribution = EXPONENTIAL;

    } else if (strcasecmp(distribution_string, "PLUMMER") == 0) {
        *distribution = PLUMMER;
        
    } else if (strcasecmp(distribution_string, "PLUMMER_SYMMETRIC") == 0) {
        *distribution = PLUMMER_SYMMETRIC;
        
    } else if (strcasecmp(distribution_string, "SLAB_1") == 0) {
        *distribution = SLAB_1;
        
    } else if (strcasecmp(distribution_string, "SLAB_2") == 0) {
        *distribution = SLAB_2;
        
    } else if (strcasecmp(distribution_string, "SPHERICAL_SHELL") == 0) {
        *distribution = SPHERICAL_SHELL;
        
    } else {
        if (rank == 0) {
            printf("[random cube example] ERROR! Undefined distribution token \"%s\". Exiting.\n",
                   distribution_string);
        }
        exit(1);
    }


    if (strcasecmp(partition_string, "RCB") == 0) {
        *partition = RCB;
        
    } else if (strcasecmp(partition_string, "HSFC") == 0) {
        *partition = HSFC;

    } else {
        if (rank == 0) {
            printf("[random cube example] ERROR! Undefined distribution token \"%s\". Exiting.\n",
                   distribution_string);
        }
        exit(1);
    }
    

    RunParams_Setup(run_params,
                    kernel, num_kernel_params, kernel_params,
                    approximation, singularity, compute_type,
                    theta, interp_degree,
                    max_per_source_leaf, max_per_target_leaf,
                    size_check_factor, beta, verbosity);

    return;
}



/*----------------------------------------------------------------------------*/
double Point_Set_Init(DISTRIBUTION distribution)
{
    if (distribution == UNIFORM) {
        return (double)random()/(double)(RAND_MAX);
        
    } else if (distribution == GAUSSIAN) {
        
        double u = (double)random()/(1.+ (double)(RAND_MAX));
        double x = 1. / sqrt(6.) * erfinv(2. * u - 1.);
	
        return x;
        
    } else if (distribution == EXPONENTIAL) {
        
        double u = (double)random()/(1.+ (double)(RAND_MAX));
        double x = -log(1. - u) / sqrt(12.);
        
        return x;

    } else {

        printf("[random cube example] ERROR! Distribution %d undefined in this "
               "context.  Exiting.\n", distribution);
        exit(1);
    }
}


/*----------------------------------------------------------------------------*/
double Point_Set(DISTRIBUTION distribution, double xmin, double xmax)
{
    double cdf_min, cdf_max;
    
    if (distribution == UNIFORM) {
        return (double)random()/(double)(RAND_MAX) * (xmax - xmin) + xmin;
        
    } else if (distribution == GAUSSIAN) {
        
        cdf_min = 0.5 * (1. + erf((xmin) * sqrt(6.)));
        cdf_max = 0.5 * (1. + erf((xmax) * sqrt(6.)));
        
        double u = (double)random()/(double)(RAND_MAX) * (cdf_max - cdf_min) + cdf_min;
        
        return 0.5 + 1. / sqrt(6.) * erfinv(2. * u - 1.);
        
    } else if (distribution == EXPONENTIAL) {
        
        cdf_min = 1 - exp(-sqrt(12) * xmin);
        if (xmax > 1) {
            cdf_max = 1;
        } else {
            cdf_max = 1 - exp(-sqrt(12) * xmax);
        }
        
        double u = (double)random()/(1. + (double)(RAND_MAX)) * (cdf_max - cdf_min) + cdf_min;
                
        return -log(1. - u) / sqrt(12.);
        
    } else {

        printf("[random cube example] ERROR! Distribution %d undefined in this "
               "context.  Exiting.\n", distribution);
        exit(1);
    }
}


/*----------------------------------------------------------------------------*/
void Point_Plummer(double R, double *x, double *y, double *z)
{
    do {
        double u = (double)random()/(1.+ (double)(RAND_MAX));
        double radius = R / sqrt(pow(u, (-2.0/3.0)) - 1.0);

        u = (double)random()/(1.+ (double)(RAND_MAX));
        double theta = acos(-1 + u * 2.0);
        
        u = (double)random()/(1.+ (double)(RAND_MAX));
        double phi = u * 2.0 * M_PI;

        *x = radius * sin(theta) * cos(phi);
        *y = radius * sin(theta) * sin(phi);
        *z = radius * cos(theta);
    } while (fabs(*x) > 100 || fabs(*y) > 100 || fabs(*z) > 100);

    return;
}


/*----------------------------------------------------------------------------*/
void Point_Plummer_Octant(double R, double *x, double *y, double *z)
{
    double u = (double)random()/(1.+ (double)(RAND_MAX));
    double radius = R / sqrt(pow(u, (-2.0/3.0)) - 1.0);

    u = (double)random()/(1.+ (double)(RAND_MAX));
    double theta = acos(u);
    
    u = (double)random()/(1.+ (double)(RAND_MAX));
    double phi = u * M_PI / 2.0;

    *x = radius * sin(theta) * cos(phi);
    *y = radius * sin(theta) * sin(phi);
    *z = radius * cos(theta);

    return;
}


/*----------------------------------------------------------------------------*/
void Point_Gaussian(double *x, double *y, double *z)
{
    double u = (double)random()/(1.+ (double)(RAND_MAX));
    *x = 1. / sqrt(6.) * erfinv(2. * u - 1.);

    u = (double)random()/(1.+ (double)(RAND_MAX));
    *y = 1. / sqrt(6.) * erfinv(2. * u - 1.);
    
    u = (double)random()/(1.+ (double)(RAND_MAX));
    *z = 1. / sqrt(6.) * erfinv(2. * u - 1.);

    return;
}


/*----------------------------------------------------------------------------*/
void Point_Exponential(double *x, double *y, double *z)
{
    double u = (double)random()/(1.+ (double)(RAND_MAX));
    *x = -log(1. - u) / sqrt(12.);

    u = (double)random()/(1.+ (double)(RAND_MAX));
    *y = -log(1. - u) / sqrt(12.);
    
    u = (double)random()/(1.+ (double)(RAND_MAX));
    *z = -log(1. - u) / sqrt(12.);

    return;
}


/*----------------------------------------------------------------------------*/
void Point_Spherical_Shell(double R, double *x, double *y, double *z)
{
    double u = 2. * M_PI * (double)random()/(1.+ (double)(RAND_MAX));
    double v = 2. * (double)random()/(1.+ (double)(RAND_MAX)) - 1.;

    *x = sqrt(1. - v * v) * cos(u);
    *y = sqrt(1. - v * v) * sin(u);
    *z = v;

    *x *= R;
    *y *= R;
    *z *= R;

    return;
}



/*----------------------------------------------------------------------------*/
void Timing_Calculate(double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
                      double time_run[4], double time_tree[13], double time_direct[4])
{
    MPI_Reduce(time_run, &time_run_glob[0], 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_run, &time_run_glob[1], 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_run, &time_run_glob[2], 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(time_direct, &time_direct_glob[0], 4, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_direct, &time_direct_glob[1], 4, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_direct, &time_direct_glob[2], 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    MPI_Reduce(time_tree, &time_tree_glob[0], 13, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[1], 13, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(time_tree, &time_tree_glob[2], 13, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    return;
}


/*----------------------------------------------------------------------------*/
void Timing_Print(double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
                  int run_direct, struct RunParams *run_params)
{
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    
    if (rank == 0) {
    
    double min_percent = 100. / time_run_glob[0][3];
    double max_percent = 100. / time_run_glob[1][3];
    double avg_percent = 100. / time_run_glob[2][3];

    double min_percent_direct = 100. / time_run_glob[0][1];
    double max_percent_direct = 100. / time_run_glob[1][1];
    double avg_percent_direct = 100. / time_run_glob[2][1];

    double min_percent_tree = 100. / time_run_glob[0][2];
    double max_percent_tree = 100. / time_run_glob[1][2];
    double avg_percent_tree = 100. / time_run_glob[2][2];

    /* Printing direct and treecode time calculations: */
    printf("[random cube example]\n");
    printf("[random cube example] ");
    printf("Treecode timing summary (all times in seconds)...\n");
    printf("[random cube example]\n");
    printf("[random cube example] ");
    printf("                                       Max                           Avg                          Max/Min\n");
    printf("[random cube example] ");
    printf("|    Total time......................  %9.3e s    (100.00%%)      %9.3e s    (100.00%%)    %8.3f \n",
                 time_run_glob[1][3], time_run_glob[2][3]/numProcs, time_run_glob[1][3]/time_run_glob[0][3]);
    printf("[random cube example] ");
    printf("|    |\n");
    printf("[random cube example] ");
    printf("|    |....Pre-process................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_run_glob[1][0],          time_run_glob[1][0] * max_percent,
                 time_run_glob[2][0]/numProcs, time_run_glob[2][0] * avg_percent,
                 time_run_glob[1][0]/time_run_glob[0][0]);

    if (run_direct == 1) {
    printf("[random cube example] ");
    printf("|    |....Directdriver...............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_run_glob[1][1],          time_run_glob[1][1] * max_percent,
                 time_run_glob[2][1]/numProcs, time_run_glob[2][1] * avg_percent,
                 time_run_glob[1][1]/time_run_glob[0][1]);
    }
    printf("[random cube example] ");
    printf("|    |....Treedriver.................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_run_glob[1][2],          time_run_glob[1][2] * max_percent,
                 time_run_glob[2][2]/numProcs, time_run_glob[2][2] * avg_percent,
                 time_run_glob[1][2]/time_run_glob[0][2]);
    printf("[random cube example]\n");
    printf("[random cube example]\n");


    if (run_direct == 1) {
    printf("[random cube example] ");
    printf("|    Directdriver....................  %9.3e s    (100.00%%)      %9.3e s    (100.00%%)    %8.3f \n",
                 time_run_glob[1][1], time_run_glob[2][1]/numProcs, time_run_glob[1][1]/time_run_glob[0][1]);

    printf("[random cube example] ");
    printf("|    |\n");

    if (numProcs > 1) {
    printf("[random cube example] ");
    printf("|    |....Communicate................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_direct_glob[1][0],          time_direct_glob[1][0] * max_percent_direct,
                 time_direct_glob[2][0]/numProcs, time_direct_glob[2][0] * avg_percent_direct,
                 time_direct_glob[1][0]/time_direct_glob[0][0]);
    }

    printf("[random cube example] ");
    printf("|    |....Compute local..............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_direct_glob[1][2],          time_direct_glob[1][2] * max_percent_direct,
                 time_direct_glob[2][2]/numProcs, time_direct_glob[2][2] * avg_percent_direct,
                 time_direct_glob[1][2]/time_direct_glob[0][2]);

    if (numProcs > 1) {
    printf("[random cube example] ");
    printf("|    |....Compute remote.............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_direct_glob[1][1],          time_direct_glob[1][1] * max_percent_direct,
                 time_direct_glob[2][1]/numProcs, time_direct_glob[2][1] * avg_percent_direct,
                 time_direct_glob[1][1]/time_direct_glob[0][1]);
    }

    if (run_params->singularity == SUBTRACTION) {
    printf("[random cube example] ");
    printf("|    |....Correct potential..........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_direct_glob[1][3],          time_direct_glob[1][3] * max_percent_direct,
                 time_direct_glob[2][3]/numProcs, time_direct_glob[2][3] * avg_percent_direct,
                 time_direct_glob[1][3]/time_direct_glob[0][3]);
    printf("[random cube example]\n");
    printf("[random cube example]\n");
    } else {
    printf("[random cube example]\n");
    printf("[random cube example]\n");
    }
    }

    printf("[random cube example] ");
    printf("|    Treedriver......................  %9.3e s    (100.00%%)      %9.3e s    (100.00%%)    %8.3f \n",
                 time_run_glob[1][2], time_run_glob[2][2]/numProcs, time_run_glob[1][2]/time_run_glob[0][2]);
    printf("[random cube example] ");
    printf("|    |\n");
    printf("[random cube example] ");
    printf("|    |....Build local tree...........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][0],          time_tree_glob[1][0] * max_percent_tree,
                 time_tree_glob[2][0]/numProcs, time_tree_glob[2][0] * avg_percent_tree,
                 time_tree_glob[1][0]/time_tree_glob[0][0]);
    printf("[random cube example] ");
    printf("|    |....Build local batches........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][1],          time_tree_glob[1][1] * max_percent_tree,
                 time_tree_glob[2][1]/numProcs, time_tree_glob[2][1] * avg_percent_tree,
                 time_tree_glob[1][1]/time_tree_glob[0][1]);
    printf("[random cube example] ");
    printf("|    |....Build local clusters.......  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][2],          time_tree_glob[1][2] * max_percent_tree,
                 time_tree_glob[2][2]/numProcs, time_tree_glob[2][2] * avg_percent_tree,
                 time_tree_glob[1][2]/time_tree_glob[0][2]);

    if (numProcs > 1) {
    printf("[random cube example] ");
    printf("|    |....Build LET..................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][3],          time_tree_glob[1][3] * max_percent_tree,
                 time_tree_glob[2][3]/numProcs, time_tree_glob[2][3] * avg_percent_tree,
                 time_tree_glob[1][3]/time_tree_glob[0][3]);
    }

    printf("[random cube example] ");
    printf("|    |....Build local lists..........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][4],          time_tree_glob[1][4] * max_percent_tree,
                 time_tree_glob[2][4]/numProcs, time_tree_glob[2][4] * avg_percent_tree,
                 time_tree_glob[1][4]/time_tree_glob[0][4]);
    printf("[random cube example] ");
    printf("|    |....Compute local..............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][5],          time_tree_glob[1][5] * max_percent_tree,
                 time_tree_glob[2][5]/numProcs, time_tree_glob[2][5] * avg_percent_tree,
                 time_tree_glob[1][5]/time_tree_glob[0][5]);

    if (numProcs > 1) {
    printf("[random cube example] ");
    printf("|    |....Build remote lists.........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][6],          time_tree_glob[1][6] * max_percent_tree,
                 time_tree_glob[2][6]/numProcs, time_tree_glob[2][6] * avg_percent_tree,
                 time_tree_glob[1][6]/time_tree_glob[0][6]);
    printf("[random cube example] ");
    printf("|    |....Compute remote.............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][7],          time_tree_glob[1][7] * max_percent_tree,
                 time_tree_glob[2][7]/numProcs, time_tree_glob[2][7] * avg_percent_tree,
                 time_tree_glob[1][7]/time_tree_glob[0][7]);
    }

    if (run_params->compute_type == CLUSTER_PARTICLE || run_params->compute_type == CLUSTER_CLUSTER) {
    printf("[random cube example] ");
    printf("|    |....Compute cp2................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][8],          time_tree_glob[1][8] * max_percent_tree,
                 time_tree_glob[2][8]/numProcs, time_tree_glob[2][8] * avg_percent_tree,
                 time_tree_glob[1][8]/time_tree_glob[0][8]);
    }

    printf("[random cube example] ");
    printf("|    |....Correct potential..........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][9],          time_tree_glob[1][9] * max_percent_tree,
                 time_tree_glob[2][9]/numProcs, time_tree_glob[2][9] * avg_percent_tree,
                 time_tree_glob[1][9]/time_tree_glob[0][9]);
    printf("[random cube example] ");
    printf("|    |....Cleanup....................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][10],          time_tree_glob[1][10] * max_percent_tree,
                 time_tree_glob[2][10]/numProcs, time_tree_glob[2][10] * avg_percent_tree,
                 time_tree_glob[1][10]/time_tree_glob[0][10]);
    printf("[random cube example]\n");
    
    if (numProcs > 1) {
    printf("[random cube example] ");
    printf("((   |....Total setup................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f ))\n",
                 time_tree_glob[1][11],          time_tree_glob[1][11] * max_percent_tree,
                 time_tree_glob[2][11]/numProcs, time_tree_glob[2][11] * avg_percent_tree,
                 time_tree_glob[1][11]/time_tree_glob[0][11]);
    printf("[random cube example] ");
    printf("((   |....Build local clusters.......  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f ))\n",
                 time_tree_glob[1][02],          time_tree_glob[1][02] * max_percent_tree,
                 time_tree_glob[2][02]/numProcs, time_tree_glob[2][02] * avg_percent_tree,
                 time_tree_glob[1][02]/time_tree_glob[0][02]);
    printf("[random cube example] ");
    printf("((   |....Total compute..............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f ))\n",
                 time_tree_glob[1][12],          time_tree_glob[1][12] * max_percent_tree,
                 time_tree_glob[2][12]/numProcs, time_tree_glob[2][12] * avg_percent_tree,
                 time_tree_glob[1][12]/time_tree_glob[0][12]);
    printf("[random cube example]\n");
    printf("[random cube example]\n");
    }
    }

    return;
}



/*----------------------------------------------------------------------------*/
void Accuracy_Calculate(double *potential_engy_glob, double *potential_engy_direct_glob,
                double *glob_inf_err, double *glob_relinf_err, double *glob_n2_err, double *glob_reln2_err,
                double *potential, double *potential_direct, int targets_num, int slice)
{
    double potential_engy = sum(potential, targets_num);
    double potential_engy_direct = sum(potential_direct, targets_num / slice);

    MPI_Reduce(&potential_engy, potential_engy_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&potential_engy_direct, potential_engy_direct_glob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    double inferr = 0.0, relinferr = 0.0, n2err = 0.0, reln2err = 0.0;
    double temp;

    for (int j = 0; j < targets_num / slice; j++) {

        temp = fabs(potential_direct[j] - potential[j*slice]);

        if (temp >= inferr) inferr = temp;

        if (fabs(potential_direct[j]) >= relinferr)
            relinferr = fabs(potential_direct[j]);

        n2err = n2err + pow(potential_direct[j] - potential[j*slice], 2.0);
        reln2err = reln2err + pow(potential_direct[j], 2.0);
    }

    MPI_Reduce(&reln2err, glob_reln2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&n2err, glob_n2_err, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&relinferr, glob_relinf_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&inferr, glob_inf_err, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    *glob_reln2_err = sqrt(fabs(*glob_n2_err / *glob_reln2_err));
    *glob_n2_err = sqrt(fabs(*glob_n2_err));
    *glob_relinf_err = *glob_inf_err / *glob_relinf_err;
    
    return;
}


/*----------------------------------------------------------------------------*/
void Accuracy_Print(double potential_engy_glob, double potential_engy_direct_glob,
                double glob_inf_err, double glob_relinf_err, double glob_n2_err, double glob_reln2_err,
                int slice)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0) {
        printf("[random cube example]                Tree potential energy:  %f\n", potential_engy_glob);
        if (slice == 1) {
        printf("[random cube example]              Direct potential energy:  %f\n", potential_engy_direct_glob);
        printf("[random cube example]\n");
        printf("[random cube example]   Absolute error for total potential:  %e\n",
               fabs(potential_engy_glob-potential_engy_direct_glob));
        printf("[random cube example]   Relative error for total potential:  %e\n",
               fabs((potential_engy_glob-potential_engy_direct_glob)/potential_engy_direct_glob));
        }
        printf("[random cube example]\n");
        printf("[random cube example] Relative inf norm error in potential:  %e \n", glob_relinf_err);
        printf("[random cube example]   Relative 2 norm error in potential:  %e \n", glob_reln2_err);
        printf("[random cube example]\n");
    }

    return;
}



/*----------------------------------------------------------------------------*/
void CSV_Print(int N, int M, struct RunParams *run_params,
               double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
               double potential_engy_glob, double potential_engy_direct_glob,
               double glob_inf_err, double glob_relinf_err, double glob_n2_err, double glob_reln2_err)
{
    int rank, numProcs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    
    if (rank == 0) {
        FILE *fp = fopen("out.csv", "a");
        fprintf(fp, "%d,%d,%d,%d,%d,%d,%d,%f,%d,%d,%d,%f,%f,"
                    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                    "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                    "%e,%e,%e,%e,%e,%e,%e,%e\n",
            N, M, numProcs, run_params->kernel, run_params->approximation, run_params->singularity,
            run_params->compute_type, run_params->theta, run_params->interp_degree,
            run_params->max_per_source_leaf, run_params->max_per_target_leaf,
            run_params->size_check_factor, run_params->beta, // 1 ends

            time_run_glob[0][0],  time_run_glob[1][0],  // min, max, avg pre-process
            time_run_glob[2][0]/numProcs,

            time_run_glob[0][1],  time_run_glob[1][1],  // min, max, avg directdriver
            time_run_glob[2][1]/numProcs,

            time_direct_glob[0][0], time_direct_glob[1][0],     // min, max, avg communicate
            time_direct_glob[2][0]/numProcs,
            time_direct_glob[0][1], time_direct_glob[1][1],     // min, max, avg compute local
            time_direct_glob[2][1]/numProcs,
            time_direct_glob[0][2], time_direct_glob[1][2],     // min, max, avg compute remote
            time_direct_glob[2][2]/numProcs,
            time_direct_glob[0][3], time_direct_glob[1][3],     // min, max, avg correct potential
            time_direct_glob[2][3]/numProcs, // 2 ends

            time_run_glob[0][2],  time_run_glob[1][2],  // min, max, avg treedriver
            time_run_glob[2][2]/numProcs,

            time_tree_glob[0][0], time_tree_glob[1][0],     // min, max, avg build local tree
            time_tree_glob[2][0]/numProcs,
            time_tree_glob[0][1], time_tree_glob[1][1],     // min, max, avg build local batches
            time_tree_glob[2][1]/numProcs,
            time_tree_glob[0][2], time_tree_glob[1][2],     // min, max, avg fill local clusters
            time_tree_glob[2][2]/numProcs,
            time_tree_glob[0][3], time_tree_glob[1][3],     // min, max, avg build LET
            time_tree_glob[2][3]/numProcs,
            time_tree_glob[0][4], time_tree_glob[1][4],     // min, max, avg build local lists
            time_tree_glob[2][4]/numProcs, // 3 ends

            time_tree_glob[0][5], time_tree_glob[1][5],     // min, max, avg compute local interations
            time_tree_glob[2][5]/numProcs,
            time_tree_glob[0][6], time_tree_glob[1][6],     // min, max, avg build remote lists
            time_tree_glob[2][6]/numProcs,
            time_tree_glob[0][7], time_tree_glob[1][7],     // min, max, avg compute remote interactions
            time_tree_glob[2][7]/numProcs,
            time_tree_glob[0][8], time_tree_glob[1][8],     // min, max, avg compute CP2 (cluster-particle)
            time_tree_glob[2][8]/numProcs,
            time_tree_glob[0][9], time_tree_glob[1][9],     // min, max, avg correct potential
            time_tree_glob[2][9]/numProcs,
            time_tree_glob[0][10], time_tree_glob[1][10],     // min, max, avg cleanup
            time_tree_glob[2][10]/numProcs,

            time_tree_glob[0][11], time_tree_glob[1][11],   // min, max, avg total setup
            time_tree_glob[2][11]/numProcs,
            time_tree_glob[0][12], time_tree_glob[1][12],   // min, max, avg total cleanup
            time_tree_glob[2][12]/numProcs,

            time_run_glob[0][3],  time_run_glob[1][3],  // min, max, avg total time
            time_run_glob[2][3]/numProcs, // 4 ends

            potential_engy_direct_glob, potential_engy_glob,
            fabs(potential_engy_direct_glob - potential_engy_glob),
            fabs((potential_engy_direct_glob - potential_engy_glob) / potential_engy_direct_glob),
            glob_inf_err, glob_relinf_err, glob_n2_err, glob_reln2_err); // 5 ends
        fclose(fp);
    }

    return;
}



/*----------------------------------------------------------------------------*/
#define erfinv_a3 -0.140543331
#define erfinv_a2 0.914624893
#define erfinv_a1 -1.645349621
#define erfinv_a0 0.886226899

#define erfinv_b4 0.012229801
#define erfinv_b3 -0.329097515
#define erfinv_b2 1.442710462
#define erfinv_b1 -2.118377725
#define erfinv_b0 1

#define erfinv_c3 1.641345311
#define erfinv_c2 3.429567803
#define erfinv_c1 -1.62490649
#define erfinv_c0 -1.970840454

#define erfinv_d2 1.637067800
#define erfinv_d1 3.543889200
#define erfinv_d0 1

double erfinv (double x)
{
    double x2, r, y;
    int  sign_x;

    if (x < -1 || x > 1)
    return NAN;

    if (x == 0)
    return 0;

    if (x > 0)
    sign_x = 1;
    else {
    sign_x = -1;
    x = -x;
    }

    if (x <= 0.7) {

    x2 = x * x;
    r =
      x * (((erfinv_a3 * x2 + erfinv_a2) * x2 + erfinv_a1) * x2 + erfinv_a0);
    r /= (((erfinv_b4 * x2 + erfinv_b3) * x2 + erfinv_b2) * x2 +
      erfinv_b1) * x2 + erfinv_b0;
    }
    else {
    y = sqrt (-log ((1 - x) / 2));
    r = (((erfinv_c3 * y + erfinv_c2) * y + erfinv_c1) * y + erfinv_c0);
    r /= ((erfinv_d2 * y + erfinv_d1) * y + erfinv_d0);
    }

    r = r * sign_x;
    x = x * sign_x;

    r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));
    r -= (erf (r) - x) / (2 / sqrt (M_PI) * exp (-r * r));
    
    return r;
}
