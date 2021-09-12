#ifndef H_SUPPORT_FUNCTIONS_H
#define H_SUPPORT_FUNCTIONS_H

#include <stdlib.h>

#include "../src/run_params/struct_run_params.h"


void Params_Parse_Readin(FILE *fp, struct RunParams **run_params, int *N, char *file_pqr,
                         int *run_direct, int *slice, double *xyz_limits, int *grid_dim);


void Timing_Calculate(double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
                double time_run[4], double time_tree[13], double time_direct[4]);
                      
void Timing_Print(double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
                int run_direct, struct RunParams *run_params);
                  
                  
void Accuracy_Calculate(double *potential_engy_glob, double *potential_engy_direct_glob,
                double *glob_inf_err, double *glob_relinf_err, double *glob_n2_err, double *glob_reln2_err,
                double *potential, double *potential_direct, int *grid_dim, int *slice);
                
void Accuracy_Print(double potential_engy_glob, double potential_engy_direct_glob,
                double glob_inf_err, double glob_relinf_err, double glob_n2_err, double glob_reln2_err,
                int *slice);
                
                
void CSV_Print(int N, int M, struct RunParams *run_params,
                double time_run_glob[3][4], double time_tree_glob[3][13], double time_direct_glob[3][4],
                double potential_engy_glob, double potential_engy_direct_glob,
                double glob_inf_err, double glob_relinf_err, double glob_n2_err, double glob_reln2_err);


#endif /* H_SUPPORT_FUNCTIONS_H */
