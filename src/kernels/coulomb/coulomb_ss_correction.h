/* Interaction Kernels */
#ifndef H_K_COULOMB_SS_CORRECTION_H
#define H_K_COULOMB_SS_CORRECTION_H
 
#include "../../run_params/struct_run_params.h"


void K_Coulomb_SS_Correction(double *potential, double *target_q,
        int numTargets, struct RunParams *run_params);


#endif /* H_K_COULOMB_SS_CORRECTION_H */
