#ifndef H_ADVANCED_TIMINGS_H
#define H_ADVANCED_TIMINGS_H

#include <stdlib.h>

#include "../run_params/struct_run_params.h"


void Timing_Calculate(double time_tree_glob[3][13], double time_tree[13], double total_time_glob[1], double total_time[1]);
                      
void Timing_Print(double time_tree_glob[3][13], double total_time_glob[1], struct RunParams *run_params);
                  
                  
#endif /* H_ADVANCED_TIMINGS_H */
