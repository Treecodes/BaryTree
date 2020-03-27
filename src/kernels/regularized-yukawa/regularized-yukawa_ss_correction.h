/* Interaction Kernels */
#ifndef H_REGULARIZED_YUKAWA_SS_CORRECTION_H
#define H_REGULARIZED_YUKAWA_SS_CORRECTION_H
 
#include "../../run_params/struct_run_params.h"


void K_RegularizedYukawa_SS_Correction(double *potential, double *target_q,
        int numTargets, struct RunParams *run_params);


#endif /* H_K_REGULARIZED_YUKAWA_SS_CORRECTION_H */
