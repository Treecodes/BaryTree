#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "../src/utilities/tools.h"
#include "../src/run_params/run_params.h"
#include "../src/run_params/struct_run_params.h"

#include "support_fns.h"


/*----------------------------------------------------------------------------*/
void Params_Parse_Readin(FILE *fp, struct RunParams **run_params, int *N, char *file_pqr,
                         int *run_direct, int *slice, double *xyz_limits, int *grid_dim)
{
    /* BaryTree params */
    int verbosity = 0;
    int interp_order = 5; 
    double theta = 0.5; 
    int max_per_source_leaf = 500;
    int max_per_target_leaf = 500; 
    double size_check_factor = 1.0;
    
    char kernel_string[256]        = "COULOMB";
    char approximation_string[256] = "LAGRANGE";
    char compute_type_string[256]  = "CLUSTER_PARTICLE";
    char run_direct_string[256]    = "OFF";

    KERNEL kernel;
    APPROXIMATION approximation;
    COMPUTE_TYPE compute_type;

    int num_kernel_params = 0;
    double kernel_params[32];

    /* random_cube_example params */
    *N = 10000; 
    slice[0] = 1;
    slice[1] = 1;
    slice[2] = 1;
    
    xyz_limits[0] = -1.;
    xyz_limits[1] =  1.;
    xyz_limits[2] = -1.;
    xyz_limits[3] =  1.;
    xyz_limits[4] = -1.;
    xyz_limits[5] =  1.;

    grid_dim[0] = 100;
    grid_dim[1] = 100;
    grid_dim[2] = 100;


    char c[256], c1[256], c2[256];

    while (fgets(c, 256, fp) != NULL) {
        sscanf(c, "%s %s", c1, c2);
    
        /* Parameters for the RunParam struct */
        if (strcmp(c1, "order") == 0) {
            interp_order = atoi(c2);

        } else if (strcmp(c1, "theta") == 0) {
            theta = atof(c2);

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
        
        } else if (strcmp(c1, "approximation") == 0) {
            strcpy(approximation_string, c2);
        
        } else if (strcmp(c1, "compute_type") == 0) {
            strcpy(compute_type_string, c2);
        
        } else if (strcmp(c1, "size_check") == 0) {
            size_check_factor = atof(c2);
        
        } else if (strcmp(c1, "verbosity") == 0) {
            verbosity = atoi(c2);

        /* Other run parameters */
        } else if (strcmp(c1, "num_sources") == 0) {
            *N = atoi(c2);

        } else if (strcmp(c1, "run_direct") == 0) {
            strcpy(run_direct_string, c2);

        } else if (strcmp(c1, "slice") == 0) {
            char *slice_string = strtok(c2, " ,");
            int num_slice = 0;
            while (slice_string != NULL) {
                slice[num_slice] = atoi(slice_string);
                num_slice++;
                slice_string = strtok(NULL, " ,");
            }
            
        } else if (strcmp(c1, "file_pqr") == 0) {
            strcpy(file_pqr, c2);

        } else if (strcmp(c1, "xyz_limits") == 0) {
            char *xyz_limits_string = strtok(c2, " ,");
            int num_xyz_limits = 0;
            while (xyz_limits_string != NULL) {
                xyz_limits[num_xyz_limits] = atof(xyz_limits_string);
                num_xyz_limits++;
                xyz_limits_string = strtok(NULL, " ,");
            }

        } else if (strcmp(c1, "grid_dim") == 0) {
            char *grid_dim_string = strtok(c2, " ,");
            int num_grid_dim = 0;
            while (grid_dim_string != NULL) {
                grid_dim[num_grid_dim] = atoi(grid_dim_string);
                num_grid_dim++;
                grid_dim_string = strtok(NULL, " ,");
            }

        } else {
            printf("[random cube example] ERROR! Undefined token \"%s\". Exiting.\n", c1);
            exit(1);
        }
    }
    

    /* Validating tokens for RunParam struct */
    if (strcasecmp(kernel_string, "COULOMB") == 0) {
        kernel = COULOMB;

    } else if (strcasecmp(kernel_string, "YUKAWA") == 0) {
        kernel = YUKAWA;
    
    } else if (strcasecmp(kernel_string, "TCF") == 0) {
        kernel = TCF;
    
    } else if (strcasecmp(kernel_string, "DCF") == 0) {
        kernel = DCF;

    } else {
        printf("[random cube example] ERROR! Undefined kernel token \"%s\". Exiting.\n",
               kernel_string);
        exit(1);
    }


    if (strcasecmp(approximation_string, "LAGRANGE") == 0) {
        approximation = LAGRANGE;

    } else if (strcasecmp(approximation_string, "HERMITE") == 0) {
        approximation = HERMITE;
    
    } else {
        printf("[random cube example] ERROR! Undefined approximation token \"%s\". Exiting.\n",
               approximation_string);
        exit(1);
    }


    if ((strcasecmp(compute_type_string, "CLUSTER_PARTICLE") == 0)
     || (strcasecmp(compute_type_string, "CLUSTER-PARTICLE") == 0)) {
        compute_type = CLUSTER_PARTICLE;

    } else if ((strcasecmp(compute_type_string, "CLUSTER_CLUSTER") == 0)
            || (strcasecmp(compute_type_string, "CLUSTER-CLUSTER") == 0)) {
        compute_type = CLUSTER_CLUSTER;

    } else {
        printf("[random cube example] ERROR! Undefined compute_type token \"%s\". Exiting.\n",
               compute_type_string);
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
        printf("[random cube example] ERROR! Undefined run direct token \"%s\". Exiting.\n",
               run_direct_string);
        exit(1);
    }
    
    

    RunParams_Setup(run_params,
                    kernel, num_kernel_params, kernel_params,
                    approximation, compute_type,
                    theta, size_check_factor, interp_order, 
                    max_per_source_leaf, max_per_target_leaf,
                    verbosity);

    return;
}


/*----------------------------------------------------------------------------*/
void Timing_Calculate(double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
                      double time_run[4], double time_tree[13], double time_direct[4])
{
    for (int i = 0; i < 4; ++i) {
        time_run_glob[0][i] = time_run[i];
        time_run_glob[1][i] = time_run[i];
        time_run_glob[2][i] = time_run[i];
        time_direct_glob[0][i] = time_direct[i];
        time_direct_glob[1][i] = time_direct[i];
        time_direct_glob[2][i] = time_direct[i];
    }

    for (int i = 0; i < 13; ++i) {
        time_tree_glob[0][i] = time_tree[i];
        time_tree_glob[1][i] = time_tree[i];
        time_tree_glob[2][i] = time_tree[i];
    }

    return;
}


/*----------------------------------------------------------------------------*/
void Timing_Print(double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
                  int run_direct, struct RunParams *run_params)
{
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
                 time_run_glob[1][3], time_run_glob[2][3], time_run_glob[1][3]/time_run_glob[0][3]);
    printf("[random cube example] ");
    printf("|    |\n");
    printf("[random cube example] ");
    printf("|    |....Pre-process................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_run_glob[1][0],          time_run_glob[1][0] * max_percent,
                 time_run_glob[2][0], time_run_glob[2][0] * avg_percent,
                 time_run_glob[1][0]/time_run_glob[0][0]);

    if (run_direct == 1) {
    printf("[random cube example] ");
    printf("|    |....Directdriver...............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_run_glob[1][1],          time_run_glob[1][1] * max_percent,
                 time_run_glob[2][1], time_run_glob[2][1] * avg_percent,
                 time_run_glob[1][1]/time_run_glob[0][1]);
    }
    printf("[random cube example] ");
    printf("|    |....Treedriver.................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_run_glob[1][2],          time_run_glob[1][2] * max_percent,
                 time_run_glob[2][2], time_run_glob[2][2] * avg_percent,
                 time_run_glob[1][2]/time_run_glob[0][2]);
    printf("[random cube example]\n");
    printf("[random cube example]\n");


    if (run_direct == 1) {
    printf("[random cube example] ");
    printf("|    Directdriver....................  %9.3e s    (100.00%%)      %9.3e s    (100.00%%)    %8.3f \n",
                 time_run_glob[1][1], time_run_glob[2][1], time_run_glob[1][1]/time_run_glob[0][1]);

    printf("[random cube example] ");
    printf("|    |\n");

    printf("[random cube example] ");
    printf("|    |....Compute local..............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_direct_glob[1][2],          time_direct_glob[1][2] * max_percent_direct,
                 time_direct_glob[2][2], time_direct_glob[2][2] * avg_percent_direct,
                 time_direct_glob[1][2]/time_direct_glob[0][2]);

    printf("[random cube example]\n");
    printf("[random cube example]\n");
    }

    printf("[random cube example] ");
    printf("|    Treedriver......................  %9.3e s    (100.00%%)      %9.3e s    (100.00%%)    %8.3f \n",
                 time_run_glob[1][2], time_run_glob[2][2], time_run_glob[1][2]/time_run_glob[0][2]);
    printf("[random cube example] ");
    printf("|    |\n");
    printf("[random cube example] ");
    printf("|    |....Build local tree...........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][0],          time_tree_glob[1][0] * max_percent_tree,
                 time_tree_glob[2][0], time_tree_glob[2][0] * avg_percent_tree,
                 time_tree_glob[1][0]/time_tree_glob[0][0]);
    printf("[random cube example] ");
    printf("|    |....Build local batches........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][1],          time_tree_glob[1][1] * max_percent_tree,
                 time_tree_glob[2][1], time_tree_glob[2][1] * avg_percent_tree,
                 time_tree_glob[1][1]/time_tree_glob[0][1]);
    printf("[random cube example] ");
    printf("|    |....Build local clusters.......  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][2],          time_tree_glob[1][2] * max_percent_tree,
                 time_tree_glob[2][2], time_tree_glob[2][2] * avg_percent_tree,
                 time_tree_glob[1][2]/time_tree_glob[0][2]);

    printf("[random cube example] ");
    printf("|    |....Build local lists..........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][4],          time_tree_glob[1][4] * max_percent_tree,
                 time_tree_glob[2][4], time_tree_glob[2][4] * avg_percent_tree,
                 time_tree_glob[1][4]/time_tree_glob[0][4]);
    printf("[random cube example] ");
    printf("|    |....Compute local..............  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][5],          time_tree_glob[1][5] * max_percent_tree,
                 time_tree_glob[2][5], time_tree_glob[2][5] * avg_percent_tree,
                 time_tree_glob[1][5]/time_tree_glob[0][5]);

    printf("[random cube example] ");
    printf("|    |....Compute cp2................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][8],          time_tree_glob[1][8] * max_percent_tree,
                 time_tree_glob[2][8], time_tree_glob[2][8] * avg_percent_tree,
                 time_tree_glob[1][8]/time_tree_glob[0][8]);

    printf("[random cube example] ");
    printf("|    |....Correct potential..........  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][9],          time_tree_glob[1][9] * max_percent_tree,
                 time_tree_glob[2][9], time_tree_glob[2][9] * avg_percent_tree,
                 time_tree_glob[1][9]/time_tree_glob[0][9]);
    printf("[random cube example] ");
    printf("|    |....Cleanup....................  %9.3e s    (%6.2f%%)      %9.3e s    (%6.2f%%)    %8.3f \n",
                 time_tree_glob[1][10],          time_tree_glob[1][10] * max_percent_tree,
                 time_tree_glob[2][10], time_tree_glob[2][10] * avg_percent_tree,
                 time_tree_glob[1][10]/time_tree_glob[0][10]);
    printf("[random cube example]\n");
    
    return;
}



/*----------------------------------------------------------------------------*/
void Accuracy_Calculate(double *potential_engy_glob, double *potential_engy_direct_glob,
                double *glob_inf_err, double *glob_relinf_err, double *glob_n2_err, double *glob_reln2_err,
                double *potential, double *potential_direct, int *grid_dim, int *slice)
{
    int targets_num = grid_dim[0]*grid_dim[1]*grid_dim[2];
    int targets_sample_num = (grid_dim[0]/slice[0]) * (grid_dim[1]/slice[1]) * (grid_dim[2]/slice[2]);

    double potential_engy = sum(potential, targets_num);
    double potential_engy_direct = sum(potential_direct, targets_sample_num);

    *potential_engy_glob = potential_engy;
    *potential_engy_direct_glob = potential_engy_direct;

    double inferr = 0.0, relinferr = 0.0, n2err = 0.0, reln2err = 0.0;
    double temp;

    for (int ix = 0; ix < grid_dim[0]/slice[0]; ix++) {                                                       
        for (int iy = 0; iy < grid_dim[1]/slice[1]; iy++) {                                                   
            for (int iz = 0; iz < grid_dim[2]/slice[2]; iz++) {                                               
                                                                                                                           
                int ii_sample = (ix * grid_dim[1]/slice[1] * grid_dim[2]/slice[2]) + (iy*grid_dim[2]/slice[2]) + iz;   
                int ii = (ix*slice[0] * grid_dim[1]*grid_dim[2]) + (iy*slice[1] * grid_dim[2]) + iz*slice[2];

                temp = fabs(potential_direct[ii_sample] - potential[ii]);

                if (temp >= inferr) inferr = temp;

                if (fabs(potential_direct[ii_sample]) >= relinferr)
                    relinferr = fabs(potential_direct[ii_sample]);

                n2err = n2err + pow(potential_direct[ii_sample] - potential[ii], 2.0);
                reln2err = reln2err + pow(potential_direct[ii_sample], 2.0);
            }
        }
    }

    *glob_reln2_err = sqrt(fabs(n2err / reln2err));
    *glob_n2_err = sqrt(fabs(n2err));
    *glob_relinf_err = inferr / relinferr;
    *glob_inf_err = inferr;
    
    return;
}


/*----------------------------------------------------------------------------*/
void Accuracy_Print(double potential_engy_glob, double potential_engy_direct_glob,
                double glob_inf_err, double glob_relinf_err, double glob_n2_err, double glob_reln2_err,
                int *slice)
{
    printf("[random cube example]                Tree potential energy:  %f\n", potential_engy_glob);
    if (slice[0]*slice[1]*slice[2] == 1) {
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

    return;
}



/*----------------------------------------------------------------------------*/
void CSV_Print(int N, int M, struct RunParams *run_params,
               double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
               double potential_engy_glob, double potential_engy_direct_glob,
               double glob_inf_err, double glob_relinf_err, double glob_n2_err, double glob_reln2_err)
{
    FILE *fp = fopen("out.csv", "a");
    fprintf(fp, "%d,%d,%d,%f,%d,%d,%d,%d,"
                "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,"
                "%e,%e,%e,%e,%e,%e,%e,%e\n",
        N, M, run_params->interp_order, run_params->theta,
        run_params->max_per_source_leaf, run_params->max_per_target_leaf, run_params->kernel,
        run_params->approximation, // 1 ends

        time_run_glob[0][0],  time_run_glob[1][0],  // min, max, avg pre-process
        time_run_glob[2][0],

        time_run_glob[0][1],  time_run_glob[1][1],  // min, max, avg directdriver
        time_run_glob[2][1],

        time_direct_glob[0][0], time_direct_glob[1][0],     // min, max, avg communicate
        time_direct_glob[2][0],
        time_direct_glob[0][1], time_direct_glob[1][1],     // min, max, avg compute local
        time_direct_glob[2][1],
        time_direct_glob[0][2], time_direct_glob[1][2],     // min, max, avg compute remote
        time_direct_glob[2][2],
        time_direct_glob[0][3], time_direct_glob[1][3],     // min, max, avg correct potential
        time_direct_glob[2][3], // 2 ends

        time_run_glob[0][2],  time_run_glob[1][2],  // min, max, avg treedriver
        time_run_glob[2][2],

        time_tree_glob[0][0], time_tree_glob[1][0],     // min, max, avg build local tree
        time_tree_glob[2][0],
        time_tree_glob[0][1], time_tree_glob[1][1],     // min, max, avg build local batches
        time_tree_glob[2][1],
        time_tree_glob[0][2], time_tree_glob[1][2],     // min, max, avg fill local clusters
        time_tree_glob[2][2],
        time_tree_glob[0][3], time_tree_glob[1][3],     // min, max, avg build LET
        time_tree_glob[2][3],
        time_tree_glob[0][4], time_tree_glob[1][4],     // min, max, avg build local lists
        time_tree_glob[2][4], // 3 ends

        time_tree_glob[0][5], time_tree_glob[1][5],     // min, max, avg compute local interations
        time_tree_glob[2][5],
        time_tree_glob[0][6], time_tree_glob[1][6],     // min, max, avg build remote lists
        time_tree_glob[2][6],
        time_tree_glob[0][7], time_tree_glob[1][7],     // min, max, avg compute remote interactions
        time_tree_glob[2][7],
        time_tree_glob[0][8], time_tree_glob[1][8],     // min, max, avg compute CP2 (cluster-particle)
        time_tree_glob[2][8],
        time_tree_glob[0][9], time_tree_glob[1][9],     // min, max, avg correct potential
        time_tree_glob[2][9],
        time_tree_glob[0][10], time_tree_glob[1][10],     // min, max, avg cleanup
        time_tree_glob[2][10],

        time_tree_glob[0][11], time_tree_glob[1][11],   // min, max, avg total setup
        time_tree_glob[2][11],
        time_tree_glob[0][12], time_tree_glob[1][12],   // min, max, avg total cleanup
        time_tree_glob[2][12],

        time_run_glob[0][3],  time_run_glob[1][3],  // min, max, avg total time
        time_run_glob[2][3], // 4 ends

        potential_engy_direct_glob, potential_engy_glob,
        fabs(potential_engy_direct_glob - potential_engy_glob),
        fabs((potential_engy_direct_glob - potential_engy_glob) / potential_engy_direct_glob),
        glob_inf_err, glob_relinf_err, glob_n2_err, glob_reln2_err); // 5 ends
    fclose(fp);

    return;
}
