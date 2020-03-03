#include <stdlib.h>
#include <stdio.h>
#include <strings.h>

#include "../src/run_params/run_params.h"
#include "../src/run_params/struct_run_params.h"

#include "support_fns.h"


void Parse_Params(FILE *fp, struct RunParams **run_params, int *N, int *M, int *run_direct, int *slice)
{
    int verbosity = 0;
    int interp_order = 5; 
    double theta = 0.5; 
    int max_per_source_leaf = 500;
    int max_per_target_leaf = 500; 
    double size_check_factor = 1.0;

    char kernel_string[256] = "COULOMB";
    char singularity_string[256] = "SKIPPING";
    char approximation_string[256] = "LAGRANGE";
    char compute_type_string[256] = "PARTICLE_CLUSTER";

    KERNEL kernel;
    SINGULARITY singularity;
    APPROXIMATION approximation;
    COMPUTE_TYPE compute_type;

    int num_kernel_params = 0;
    double kernel_params[32];

    *N = 10000; 
    *M = 10000; 
    *run_direct = 0;
    *slice = 1;


    char c[256], c1[256], c2[256];

    while (fgets(c, 256, fp) != NULL) {
        sscanf(c, "%s %s", c1, c2);
    
        // Parameters for the RunParam struct
        if (strncmp(c1, "order", 5) == 0) {
            interp_order = atoi(c2);

        } else if (strncmp(c1, "theta", 5) == 0) {
            theta = atof(c2);

        } else if (strncmp(c1, "max_per_source_leaf", 19) == 0) {
            max_per_source_leaf = atoi(c2);

        } else if (strncmp(c1, "max_per_target_leaf", 19) == 0) {
            max_per_target_leaf = atoi(c2);

        } else if (strncmp(c1, "kernel_name", 11) == 0) {
            strcpy(kernel_string, c2);
        
        } else if (strncmp(c1, "kernel_params", 13) == 0) {
            char *k_param = strtok(c2, " ,");
            while (k_param != NULL) {
                kernel_params[num_kernel_params] = atof(k_param);
                num_kernel_params++;
                k_param = strtok(NULL, " ,");
            }
        
        } else if (strncmp(c1, "singularity", 11) == 0) {
            strcpy(singularity_string, c2);
        
        } else if (strncmp(c1, "approximation", 13) == 0) {
            strcpy(approximation_string, c2);
        
        } else if (strncmp(c1, "compute_type", 12) == 0) {
            strcpy(compute_type_string, c2);
        
        } else if (strncmp(c1, "size_check", 10) == 0) {
            size_check_factor = atof(c2);
        
        } else if (strncmp(c1, "verbosity", 9) == 0) {
            verbosity = atoi(c2);

        // Other run parameters
        } else if (strncmp(c1, "num_particles", 13) == 0) {
            *N = atoi(c2);
            *M = atoi(c2);

        } else if (strncmp(c1, "num_sources", 11) == 0) {
            *N = atoi(c2);

        } else if (strncmp(c1, "num_targets", 11) == 0) {
            *M = atoi(c2);

        } else if (strncmp(c1, "run_direct", 10) == 0) {
            *run_direct = atoi(c2);

        } else if (strncmp(c1, "slice", 5) == 0) {
            *slice = atoi(c2);

        } else {
            fprintf(stderr, "ERROR! Undefined token \"%s\". Exiting.\n", c1);
            exit(1);
        }
    }


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
    
    } else if (strcasecmp(kernel_string, "TCF") == 0) {
        kernel = TCF;
    
    } else if (strcasecmp(kernel_string, "DCF") == 0) {
        kernel = DCF;
    
    } else {
        fprintf(stderr, "ERROR! Undefined kernel token \"%s\". Exiting.\n", kernel_string);
        exit(1);
    }


    if (strcasecmp(singularity_string, "SKIPPING") == 0) {
        singularity = SKIPPING;

    } else if (strcasecmp(singularity_string, "SUBTRACTION") == 0) {
        singularity = SUBTRACTION;

    } else {
        fprintf(stderr, "ERROR! Undefined singularity token \"%s\". Exiting.\n", singularity_string);
        exit(1);
    }


    if (strcasecmp(approximation_string, "LAGRANGE") == 0) {
        approximation = LAGRANGE;

    } else if (strcasecmp(approximation_string, "HERMITE") == 0) {
        approximation = HERMITE;
    
    } else {
        fprintf(stderr, "ERROR! Undefined approximation token \"%s\". Exiting.\n", approximation_string);
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
        fprintf(stderr, "ERROR! Undefined compute_type token \"%s\". Exiting.\n", compute_type_string);
        exit(1);
    }


    RunParams_Setup(run_params,
                    kernel, num_kernel_params, kernel_params,
                    approximation, singularity, compute_type,
                    theta, size_check_factor, interp_order, 
                    max_per_source_leaf, max_per_target_leaf,
                    verbosity);

    return;
}
