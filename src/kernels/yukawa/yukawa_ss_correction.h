/* Interaction Kernels */
#ifndef H_K_YUKAWA_SS_CORRECTION_H
#define H_K_YUKAWA_SS_CORRECTION_H
 
#include "../../struct_kernel.h"


void K_Yukawa_SS_Correction(double *potential, double *target_q,
        int numTargets, struct kernel *kernel);


#endif /* H_K_YUKAWA_SS_CORRECTION_H */
