#ifndef H_SUPPORT_FUNCTIONS_H
#define H_SUPPORT_FUNCTIONS_H

#include <stdlib.h>

#include "../src/run_params/struct_run_params.h"


typedef enum DISTRIBUTION
{
    NO_DISTRIBUTION,
    UNIFORM,
    GAUSSIAN,
    EXPONENTIAL,
    PLUMMER,
    PLUMMER_SYMMETRIC,
    SLAB_1,
    SLAB_2,
    SPHERICAL_SHELL
} DISTRIBUTION;

typedef enum PARTITION
{
    NO_PARTITION,
    RCB,
    HSFC
} PARTITION;


void Params_Parse(FILE *fp, struct RunParams **run_params, int *N, int *M, int *run_direct, int *slice,
                double *xyz_limits, DISTRIBUTION *distribution, PARTITION *partition);


double Point_Set_Init(DISTRIBUTION distribution);

double Point_Set(DISTRIBUTION distribution, double xmin, double xmax);

void Point_Plummer(double R, double *x, double *y, double *z);

void Point_Plummer_Octant(double R, double *x, double *y, double *z);

void Point_Gaussian(double *x, double *y, double *z);

void Point_Exponential(double *x, double *y, double *z);

void Point_Spherical_Shell(double R, double *x, double *y, double *z);


void Timing_Calculate(double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
                double time_run[4], double time_tree[13], double time_direct[4]);
                      
void Timing_Print(double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
                int run_direct, struct RunParams *run_params);
                  
                  
void Accuracy_Calculate(double *potential_engy_glob, double *potential_engy_direct_glob,
                double *glob_inf_err, double *glob_relinf_err, double *glob_n2_err, double *glob_reln2_err,
                double *potential, double *potential_direct, int targets_num, int slice);
                
void Accuracy_Print(double potential_engy_glob, double potential_engy_direct_glob,
                double glob_inf_err, double glob_relinf_err, double glob_n2_err, double glob_reln2_err,
                int slice);
                
                
void CSV_Print(int N, int M, struct RunParams *run_params,
                double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
                double potential_engy_glob, double potential_engy_direct_glob,
                double glob_inf_err, double glob_relinf_err, double glob_n2_err, double glob_reln2_err);


#endif /* H_SUPPORT_FUNCTIONS_H */
